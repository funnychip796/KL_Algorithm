# Kernighan-Lin Partitioning Algorithm (KL Algorithm)

This is a C implementation of the Kernighan-Lin (KL) heuristic algorithm for graph and hypergraph bipartitioning. The goal of this algorithm is to divide a set of nodes into two equal-sized partitions while minimizing the total weight of cut edges (Cut Size) between them.

## Directory Structure

- [main.c](./main.c): The core KL algorithm program containing input parsing, iteration simulation, cut size calculation, and the main entry point.
- [testdata/](./testdata/): UCLA physical design benchmark format files (`.nodes.txt` and `.nets.txt`) used for testing.
- [run_tests.sh](./run_tests.sh): Automated build and test runner script.
- [.gitignore](./.gitignore): Specifies intentionally untracked files to ignore.

## Quick Start

After cloning the repository, you can compile the code and automatically test it using all datasets by running a single command:

```bash
# 1. Give execution permission to the test script (if needed)
chmod +x run_tests.sh

# 2. Run the test script
./run_tests.sh
```

### What the script does:
1. Automatically detects an available C compiler (checks `gcc` -> `clang` -> `cc`).
2. Compiles `main.c` to an executable `./kl_partition` with `-O3` optimization flags.
3. Scans the `testdata/` directory for matching pairs of `.nodes.txt` and `.nets.txt` benchmark files.
4. Performs partition optimization on each pair, prints initial and final cut sizes, computes reduction metrics, and shows the results in a formatted console table.

---

## Manual Compilation & Run

If you wish to compile and run a single test dataset manually, follow these steps:

### 1. Compile the code

Compile using `gcc`:
```bash
gcc main.c -o kl_partition
```

### 2. Run a test case

Pass the `.nodes.txt` and `.nets.txt` paths as arguments to the executable:
```bash
./kl_partition testdata/spp_N151_E192_R8_232.nodes.txt testdata/spp_N151_E192_R8_232.nets.txt
```

### 3. Input Data Format

The application accepts input in a simplified UCLA benchmark format:
* **`.nodes.txt`**: Specifies the total number of nodes (`NumNodes`) and defines all node names.
* **`.nets.txt`**: Defines the hyperedges, their degree (`NetDegree`), and the list of pin names connected to each net.
