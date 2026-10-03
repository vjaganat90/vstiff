(** [canon.exe file.ml] prints the implementation [file.ml] without comments,
    doc comments or its own layout, so that only a code change alters the
    output. The rules in [sync/dune] diff it against the promoted [.canon]
    files: see verif/README.md. *)
