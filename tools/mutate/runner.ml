open Outcome

type config = { root : string; workers : int; timeout_factor : float; keep : bool }
type timing = { baseline : float; timeout : float; wall : float }
type phase = Checking | Testing
type job = { mutant : Mutant.t; phase : phase; child : Child.t }
type worker = { dir : string; log : string; mutable job : job option }

(* Release profile: a warning that a mutant causes, such as an unused variable, is not an error there. No dune cache: it
   would keep the objects of every mutant, and could replay a test result instead of running the test. *)
let dune goal = [ "dune"; "build"; "--root"; "."; "--profile"; "release"; "--cache=disabled"; goal ]

(* Starting dune takes a moment too: a floor keeps a fast project from timing out on jitter. *)
let min_timeout = 5.
let ( let* ) = Result.bind

let excerpt log =
  let lines = String.split_on_char '\n' (Files.read log) in
  String.concat "\n" (List.filteri (fun i _ -> i < 20) lines)

let all_ok children = Array.for_all (fun child -> Child.poll child = Some (Unix.WEXITED 0)) children

let run config mutants say on_outcome =
  let started = Child.now () in
  let scratch = Files.temp_dir "vstiff-mutate-" in
  let interrupted = ref false in
  let on_signal = Sys.Signal_handle (fun _ -> interrupted := true) in
  let previous = List.map (fun signal -> (signal, Sys.signal signal on_signal)) [ Sys.sigint; Sys.sigterm ] in
  let workers =
    Array.init config.workers (fun i ->
        let name = Printf.sprintf "w%d" i in
        { dir = Filename.concat scratch name; log = Filename.concat scratch (name ^ ".log"); job = None })
  in
  let finally () =
    Array.iter (fun w -> Option.iter (fun j -> Child.kill j.child) w.job) workers;
    List.iter (fun (signal, behaviour) -> Sys.set_signal signal behaviour) previous;
    if config.keep then say ("the scratch copies are in " ^ scratch) else Files.remove_tree scratch
  in
  let wait children =
    while (not !interrupted) && Array.exists (fun child -> Child.poll child = None) children do
      Child.pause ()
    done;
    if !interrupted then Array.iter Child.kill children
  in
  let start w goal = Child.spawn ~dir:w.dir ~log:w.log (dune goal) in
  let in_copy w file = Filename.concat w.dir file in
  Fun.protect ~finally (fun () ->
      say (Printf.sprintf "copying %s into %d scratch copies under %s" config.root config.workers scratch);
      Array.iter (fun w -> Files.copy_tree ~src:config.root ~dst:w.dir) workers;
      (* The unmodified copies pass: otherwise every mutant would count as killed. Each copy also builds once, so that
         the mutants after it rebuild what they touch. The slowest run is the baseline. *)
      say "running the tests on the unmodified copies";
      let baselines = Array.map (fun w -> start w "@runtest") workers in
      wait baselines;
      let* () = if !interrupted then Error "interrupted" else Ok () in
      let* () =
        match Array.find_index (fun child -> Child.poll child <> Some (Unix.WEXITED 0)) baselines with
        | None -> Ok ()
        | Some i -> Error ("the tests fail on the unmodified repository:\n" ^ excerpt workers.(i).log)
      in
      let baseline = Array.fold_left (fun slowest child -> Float.max slowest (Child.elapsed child)) 0. baselines in
      let limit = Float.max min_timeout (config.timeout_factor *. baseline) in
      (* Printing a module is how a mutant is written, so the printer must not change what the code does. *)
      let files = List.sort_uniq compare (List.map (fun (m : Mutant.t) -> m.file) mutants) in
      let originals = List.map (fun file -> (file, Files.read (in_copy workers.(0) file))) files in
      let* printed =
        List.fold_right
          (fun (file, text) acc ->
            let* rest = acc in
            let* text' = Mutant.reprint ~file text in
            Ok ((file, text') :: rest))
          originals (Ok [])
      in
      say "running the tests on the modules printed unchanged";
      List.iter (fun (file, text) -> Files.write (in_copy workers.(0) file) ~contents:text) printed;
      let check = start workers.(0) "@runtest" in
      wait [| check |];
      List.iter (fun (file, text) -> Files.write (in_copy workers.(0) file) ~contents:text) originals;
      let* () = if !interrupted then Error "interrupted" else Ok () in
      let* () =
        if all_ok [| check |] then Ok ()
        else Error ("printing the modules unchanged breaks the build or the tests:\n" ^ excerpt workers.(0).log)
      in
      say (Printf.sprintf "baseline %.1f s, timeout %.1f s" baseline limit);
      let queue = ref mutants and pending = ref (List.length mutants) in
      let finish w j outcome =
        Files.write (in_copy w j.mutant.file) ~contents:(List.assoc j.mutant.file originals);
        w.job <- None;
        decr pending;
        on_outcome j.mutant outcome
      in
      let begin_job w (m : Mutant.t) =
        match Mutant.apply ~file:m.file (List.assoc m.file originals) ~id:(Mutant.id m) with
        | Error message ->
            decr pending;
            on_outcome m (Stillborn ("cannot be applied: " ^ message))
        | Ok source ->
            Files.write (in_copy w m.file) ~contents:source;
            w.job <- Some { mutant = m; phase = Checking; child = start w "@check" }
      in
      let step w =
        match w.job with
        | None -> (
            match !queue with
            | [] -> ()
            | m :: rest ->
                queue := rest;
                begin_job w m)
        | Some j -> (
            match Child.poll j.child with
            | None ->
                if Child.elapsed j.child > limit then (
                  Child.kill j.child;
                  finish w j timeout)
            | Some status -> (
                match (j.phase, status) with
                | Checking, Unix.WEXITED 0 -> w.job <- Some { j with phase = Testing; child = start w "@runtest" }
                | Checking, _ -> finish w j (stillborn ~log:(Files.read w.log))
                | Testing, Unix.WEXITED 0 -> finish w j Survived
                | Testing, _ -> finish w j (killed ~log:(Files.read w.log))))
      in
      while !pending > 0 && not !interrupted do
        Array.iter step workers;
        Child.pause ()
      done;
      if !interrupted then Error "interrupted" else Ok { baseline; timeout = limit; wall = Child.now () -. started })
