.PHONY: all build-rust

all: release

version := 0.2.4

release_targets := aarch64-apple-darwin \
x86_64-pc-windows-gnu


release:
	cargo build --bin opengm_rts --release --target aarch64-apple-darwin;
	cp target/aarch64-apple-darwin/release/opengm_rts ./dist/opengm_rts_darwin_aarch64_v${version};
	cargo build --bin opengm_rts --release --target x86_64-pc-windows-gnu;
	cp target/x86_64-pc-windows-gnu/release/opengm_rts.exe ./dist/opengm_rts_win_x86_64_v${version}.exe;
