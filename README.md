# Pidgeot

because pose graph optimization = pgo = pidgeot

## dependencies

needs both clang 17+ for compiling and Ninja for the generator

```bash
# clang
# follow https://apt.llvm.org/
# ninja
sudo apt install ninja
```

## build

```bash
make
```

that handles ensuring cxx compiler is clang and the generator is ninja. check the makefile to change versions for your setup.

im making this change so i can commit that dogleg works
