.PHONY: make clean
cpp_dir=cpp_util
cpp_build_dir=build

make:
	cd $(cpp_dir); cmake -B $(cpp_build_dir) -S .; make -C $(cpp_build_dir)

clean:
	cd $(cpp_dir); rm -rf $(cpp_build_dir)
