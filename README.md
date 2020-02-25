# PositDebug
A debugger to detect numerical errors in applications using posits. 


## How to build

1. Get llvm and clang version 9
```

  wget http://releases.llvm.org/9.0.0/llvm-9.0.0.src.tar.xz

  wget http://releases.llvm.org/9.0.0/cfe-9.0.0.src.tar.xz
```

2. Build llvm and clang

```
      tar -xvf llvm-9.0.0.src.tar.xz
      mv llvm-9.0.0.src/ llvm
      tar -xvf cfe-9.0.0.src.tar.xz
      mv cfe-9.0.0.src/* clang
      mkdir build
      cd build
      cmake -DLLVM_ENABLE_PROJECTS=clang -G "Unix Makefiles" ../llvm
      make -j8

```

3. Set env variable LLVM_HOME to the LLVM build directory
```
  export LLVM_HOME=<path to LLVM build>
```

4. Clone PositDebug git repo.
```
  git clone https://github.com/rutgers-apl/PositDebug.git

```

5. If your compiler does not support C++11 by default, add the following line to llvm-pass/PSan/CMakefile

```
  target_compile_feature(PSanitizer PRIVATE cxx_range_for cxx_auto_type)

```

otherwise, use the followng line

```
        target_compile_features(PSanitizer PRIVATE )

```

6. Build the PSan LLVM pass

```
  cd PositDebug/llvm-pass
  mkdir build
  cd build
  cmake ../
  make

```

7. Download and build the SoftPosit library with PositDebug Extensions

```
   git clone https://gitlab.com/cerlane/SoftPosit.git
   cd SoftPosit/build/Linux-x86_64-GCC/
   make

```
  Finally return back to the top-level PositDebug directory

8. Try out the tests in regression_tests directory. Set the following environment variables

```
  export LLVM_HOME=<LLVM build directory>

  export PATH=<LLVM build directory>/bin:$PATH

  export PD_HOME=<path to the PositDebug github checkout>

  export LD_LIBRARY_PATH=$PD_HOME/runtime/obj/

  export SOFTPOSIT_HOME=$PD_HOME/SoftPosit/

  export LD_LIBRARY_PATH=$PD_HOME/runtime/obj:$LD_LIBRARY_PATH
  

```

and then,
```
  make

  ./diff-root-simple.pd.o

```

It should report that there is one instance of  more than 55 bits of error in error.log
      
