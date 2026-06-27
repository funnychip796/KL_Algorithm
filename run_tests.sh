#!/bin/bash

set -e

GREEN=$(printf '\033[0;32m')
YELLOW=$(printf '\033[1;33m')
RED=$(printf '\033[0;31m')
NC=$(printf '\033[0m')

echo -e "${YELLOW}=== KL Partitioning Algorithm: Build and Test Runner ===${NC}"

if command -v gcc >/dev/null 2>&1; then
    COMPILER="gcc"
elif command -v clang >/dev/null 2>&1; then
    COMPILER="clang"
elif command -v cc >/dev/null 2>&1; then
    COMPILER="cc"
else
    echo -e "${RED}Error: No C compiler (gcc, clang, or cc) found in PATH.${NC}"
    exit 1
fi

echo -e "Using compiler: ${GREEN}${COMPILER}${NC}"

echo -e "Compiling main.c..."
$COMPILER main.c -o kl_partition

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Compilation successful! Executable created: ./kl_partition${NC}\n"
else
    echo -e "${RED}Compilation failed!${NC}"
    exit 1
fi

TESTDATA_DIR="./testdata"
if [ ! -d "$TESTDATA_DIR" ]; then
    echo -e "${RED}Error: testdata directory not found at $TESTDATA_DIR${NC}"
    exit 1
fi

echo -e "Running tests..."
printf "%-30s | %-12s | %-10s | %-15s | %-6s | %-8s\n" "Test Case" "Initial Cut" "Final Cut" "Reduction (%)" "Passes" "Status"
printf "%-30s-|-%-12s-|-%-10s-|-%-15s-|-%-6s-|-%-8s\n" "------------------------------" "------------" "----------" "---------------" "------" "--------"

for node_file in "$TESTDATA_DIR"/*.nodes.txt; do
    if [ ! -e "$node_file" ]; then
        echo -e "${RED}No .nodes.txt test files found in $TESTDATA_DIR.${NC}"
        break
    fi

    net_file="${node_file%.nodes.txt}.nets.txt"

    if [ ! -f "$net_file" ]; then
        echo -e "${YELLOW}Warning: Matching nets file not found for $(basename "$node_file")${NC}"
        continue
    fi

    case_name=$(basename "$node_file" .nodes.txt)

    output=$(./kl_partition "$node_file" "$net_file" 2>&1)

    if [ $? -ne 0 ]; then
        printf "%-30s | %-12s | %-10s | %-15s | %-6s | %-8s\n" "$case_name" "-" "-" "-" "-" "${RED}FAILED${NC}"
        continue
    fi

    initial_cut=$(echo "$output" | grep "^Initial cut:" | awk -F':' '{print $2}' | tr -d ' ')
    final_cut=$(echo "$output" | grep "^Final cut:" | awk -F':' '{print $2}' | tr -d ' ')
    passes=$(echo "$output" | grep "^Passes:" | awk -F':' '{print $2}' | tr -d ' ')

    if [ -z "$initial_cut" ] || [ -z "$final_cut" ] || [ -z "$passes" ]; then
        printf "%-30s | %-12s | %-10s | %-15s | %-6s | %-8s\n" "$case_name" "-" "-" "-" "-" "${RED}PARSE ERR${NC}"
    else
        reduction=$((initial_cut - final_cut))
        if [ "$initial_cut" -gt 0 ]; then
            pct_reduction=$(awk -v i="$initial_cut" -v r="$reduction" 'BEGIN { printf "%.1f%%", (r/i)*100 }')
        else
            pct_reduction="0.0%"
        fi
        printf "%-30s | %-12s | %-10s | %-15s | %-6s | %-8s\n" "$case_name" "$initial_cut" "$final_cut" "$reduction ($pct_reduction)" "$passes" "${GREEN}OK${NC}"
    fi
done

echo -e "\n${GREEN}All tests completed successfully!${NC}"
echo -e "To run a single test case manually:"
echo -e "  ./kl_partition <path_to_nodes_file> <path_to_nets_file>"
