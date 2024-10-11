#!/bin/bash

# automate_tracing.sh

# Step 1: Compile the program
echo "Compiling the program..."
./compile.sh

# Step 2: Run GDB with the command file
echo "Running GDB to trace function calls..."
gdb -batch -x gdb_commands.gdb ./pinocchio.x

# Step 3: Notify completion
echo "Function call tracing completed. Output saved to function_trace.log."
