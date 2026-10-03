let read path = In_channel.with_open_bin path In_channel.input_all
let exists = Sys.file_exists
let write path ~contents = Out_channel.with_open_bin path (fun channel -> Out_channel.output_string channel contents)

let append path ~contents =
  Out_channel.with_open_gen [ Open_wronly; Open_creat; Open_append; Open_binary ] 0o644 path (fun channel ->
      Out_channel.output_string channel contents)

let skipped name = name <> "" && (name.[0] = '.' || name.[0] = '_')

let rec copy_tree ~src ~dst =
  Sys.mkdir dst 0o755;
  Array.iter
    (fun name ->
      let from = Filename.concat src name and into = Filename.concat dst name in
      if skipped name then ()
      else if Sys.is_directory from then copy_tree ~src:from ~dst:into
      else write into ~contents:(read from))
    (Sys.readdir src)

let temp_dir prefix = Filename.temp_dir prefix ""

(* lstat does not follow links, so a link to a directory is unlinked and its target survives. *)
let rec remove_tree path =
  match (Unix.lstat path).st_kind with
  | Unix.S_DIR ->
      Unix.chmod path 0o700;
      Array.iter (fun name -> remove_tree (Filename.concat path name)) (Sys.readdir path);
      Unix.rmdir path
  | _ -> Sys.remove path
