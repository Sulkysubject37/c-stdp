#!/usr/bin/env bash
set -e

# Use the virtual environment python
PYTHON="casual-stdp/bin/python"

if [ ! -f "$PYTHON" ]; then
    echo "Python venv not found at $PYTHON. Falling back to python3"
    PYTHON="python3"
fi

echo "Running Unit Tests with $PYTHON..."
$PYTHON -m pytest tests/unit/

echo "Unit tests passed."
exit 0
