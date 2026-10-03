type t = { pid : int; started : float; mutable finished : (Unix.process_status * float) option }

let now = Unix.gettimeofday
let pause () = Unix.sleepf 0.05

(* fork, not create_process: the child has to call setsid before it executes. A failed exec ends the child with
   status 127, which is what a shell reports for a missing command, and _exit skips the parent's at_exit handlers. *)
let spawn ~dir ~log argv =
  flush_all ();
  let started = now () in
  match Unix.fork () with
  | 0 -> (
      try
        ignore (Unix.setsid ());
        Unix.chdir dir;
        let out = Unix.openfile log [ Unix.O_WRONLY; Unix.O_CREAT; Unix.O_TRUNC ] 0o644 in
        Unix.dup2 out Unix.stdout;
        Unix.dup2 out Unix.stderr;
        Unix.close out;
        let null = Unix.openfile "/dev/null" [ Unix.O_RDONLY ] 0 in
        Unix.dup2 null Unix.stdin;
        Unix.close null;
        Unix.execvp (List.hd argv) (Array.of_list argv)
      with _ -> Unix._exit 127)
  | pid -> { pid; started; finished = None }

let poll child =
  match child.finished with
  | Some (status, _) -> Some status
  | None -> (
      match Unix.waitpid [ Unix.WNOHANG ] child.pid with
      | 0, _ -> None
      | _, status ->
          child.finished <- Some (status, now () -. child.started);
          Some status)

let elapsed child = match child.finished with Some (_, time) -> time | None -> now () -. child.started

(* Dune puts every action in a process group of its own, so killing dune's group would orphan a test that is still
   running. SIGTERM makes dune stop what it started and exit; only if it does not within the grace period is the
   group killed. *)
let grace = 10.

let kill child =
  if child.finished = None then (
    (try Unix.kill child.pid Sys.sigterm with Unix.Unix_error _ -> ());
    let deadline = now () +. grace in
    while poll child = None && now () < deadline do
      pause ()
    done;
    if child.finished = None then (
      (try Unix.kill (-child.pid) Sys.sigkill with Unix.Unix_error _ -> ());
      (try ignore (Unix.waitpid [] child.pid) with Unix.Unix_error _ -> ());
      child.finished <- Some (Unix.WSIGNALED Sys.sigkill, now () -. child.started)))
