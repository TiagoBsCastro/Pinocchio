# Convenience top-level Makefile
#
# Common tasks:
#   make            # build the code (delegates to src)
#   make pinocchio  # explicitly build the pinocchio.x binary
#   make run_planner# build the run_planner helper
#   make clean      # clean build artifacts
#   make docs       # generate HTML API docs into ./docs/build/html
#   make docs-clean # remove generated docs

.PHONY: all pinocchio run_planner clean help docs docs-clean

# Default: build the main binary
all:
	$(MAKE) -C src

pinocchio:
	$(MAKE) -C src pinocchio

run_planner:
	$(MAKE) -C src run_planner

clean:
	$(MAKE) -C src clean

help:
	@echo "Available targets:" && \
	echo "  all         - Build the main binary (default)" && \
	echo "  pinocchio   - Build pinocchio.x explicitly" && \
	echo "  run_planner - Build the run_planner helper" && \
	echo "  clean       - Clean build artifacts" && \
	echo "  docs        - Generate HTML API documentation (output: docs/build/html)" && \
	echo "  docs-clean  - Remove generated documentation"

# Delegate to src/Makefile for docs; src runs doxygen from ./docs
docs:
	$(MAKE) -C src docs

docs-clean:
	$(MAKE) -C src docs-clean
