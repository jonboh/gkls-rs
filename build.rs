fn main() {
    // Tell Cargo that if the given file changes, to rerun this build script.
    // println!("cargo:rerun-if-changed=src/hello.c");
    // Use the `cc` crate to build a C file and statically link it.
    cc::Build::new()
        .file("gkls/C/Dp/Src/gkls.c")
        .file("gkls/C/Dp/Src/rnd_gen.c")
        .compile("gkls");

    println!("cargo:rerun-if-changed=gkls/C/Dp/Src/gkls.h");
    println!("cargo:rerun-if-changed=gkls/C/Dp/Src/gkls.c");
    println!("cargo:rerun-if-changed=gkls/C/Dp/Src/rnd_gen.h");
    println!("cargo:rerun-if-changed=gkls/C/Dp/Src/rnd_gen.c");
}
