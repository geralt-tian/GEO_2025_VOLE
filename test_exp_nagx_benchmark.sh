#!/bin/bash

# exp_nagx benchmark script - test from 2^10 to 2^18
# Measure total communication of Alice and Bob, and their respective times

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PORT=32000
BINARY="${SCRIPT_DIR}/build/bin/exp_nagx-cheetah"

echo "=========================================="
echo "exp_nagx Benchmark - VOLE Protocol"
echo "=========================================="
printf "%-17s | %-10s | %-10s | %-19s | %-17s\n" \
    "Elements" "Alice Time" "Bob Time" "Total Comm (bytes)" "Throughput (Exp/s)"
echo "------------------|------------|------------|---------------------|-------------------"

for power in {10..18}
do
    DIM=$((2**power))

    # Clean up previous logs
    rm -f /tmp/alice_${DIM}.log /tmp/bob_${DIM}.log

    # Start Alice (server)
    timeout 300 $BINARY r=1 d=$DIM nt=1 p=$PORT > /tmp/alice_${DIM}.log 2>&1 &
    ALICE_PID=$!

    # Wait for server to start
    sleep 2

    # Start Bob (client)
    timeout 300 $BINARY r=2 d=$DIM nt=1 p=$PORT ip=127.0.0.1 > /tmp/bob_${DIM}.log 2>&1
    BOB_EXIT=$?

    # Wait for Alice to complete
    wait $ALICE_PID 2>/dev/null
    ALICE_EXIT=$?

    # Extract data from logs
    ALICE_TIME=$(grep "Exp Time" /tmp/alice_${DIM}.log 2>/dev/null | awk '{print $3}')
    ALICE_COMM=$(grep "Exp Bytes Sent" /tmp/alice_${DIM}.log 2>/dev/null | awk '{print $4}')

    BOB_TIME=$(grep "Exp Time" /tmp/bob_${DIM}.log 2>/dev/null | awk '{print $3}')
    BOB_COMM=$(grep "Exp Bytes Sent" /tmp/bob_${DIM}.log 2>/dev/null | awk '{print $4}')

    # Check if data is valid
    if [ -z "$ALICE_TIME" ] || [ -z "$BOB_TIME" ] || [ -z "$ALICE_COMM" ] || [ -z "$BOB_COMM" ]; then
        printf "2^%-2d (%8d)   | FAILED - Check logs: /tmp/alice_${DIM}.log /tmp/bob_${DIM}.log\n" $power $DIM
        continue
    fi

    # Calculate total communication
    TOTAL_COMM=$((ALICE_COMM + BOB_COMM))

    # Calculate throughput (using the larger time)
    if [ "$ALICE_TIME" -gt "$BOB_TIME" ] 2>/dev/null; then
        MAX_TIME=$ALICE_TIME
    else
        MAX_TIME=$BOB_TIME
    fi

    THROUGHPUT=$(echo "scale=2; $DIM * 1000 / $MAX_TIME" | bc)

    # Output results
    printf "2^%-2d (%8d)   | %7s ms | %7s ms | %19d | %17.0f\n" \
        $power $DIM "$ALICE_TIME" "$BOB_TIME" $TOTAL_COMM $THROUGHPUT

    # Clean up processes
    pkill -P $ALICE_PID 2>/dev/null

    # Brief pause to avoid port conflicts
    sleep 1
done

echo "=========================================="
echo "Benchmark completed!"
echo "Detailed logs: /tmp/alice_*.log and /tmp/bob_*.log"
