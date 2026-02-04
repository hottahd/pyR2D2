.PHONY: make clean
fortran_dir=fortran_util/fortran_src
cpp_dir=cpp_util
cpp_build_dir = build
make:
	cd $(fortran_dir); make
	cd $(cpp_dir); cmake -B $(cpp_build_dir) -S .; make -C $(cpp_build_dir)	
clean:
	cd $(fortran_dir); make clean
	cd $(cpp_dir); rm -rf $(cpp_build_dir)