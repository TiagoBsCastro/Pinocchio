# gdb_commands.gdb

# Load the Python script
source trace_calls.py

# Enable logging
set logging file function_trace.log
set logging redirect on
set logging on

# Enable function call tracing
tracecalls

# Run the program
run parameter_file

# Exit GDB after execution
quit

