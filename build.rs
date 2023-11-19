fn main() {
    if cfg!(feature = "test_cbinding") {
        cc::Build::new()
            .file("gkls/C/Dp/Src/gkls.c")
            .file("gkls/C/Dp/Src/rnd_gen.c")
            .compile("gkls");

        println!("cargo:rerun-if-changed=gkls/C/Dp/Src/gkls.h");
        println!("cargo:rerun-if-changed=gkls/C/Dp/Src/gkls.c");
        println!("cargo:rerun-if-changed=gkls/C/Dp/Src/rnd_gen.h");
        println!("cargo:rerun-if-changed=gkls/C/Dp/Src/rnd_gen.c");
    }
}
