#!/bin/bash

# exp_nagx benchmark script - test from 2^10 to 2^18
# 统计 Alice 和 Bob 的通信量总和，以及各自的时间

PORT=32000
BINARY="/root/SEAF-VOLE-AF/build/bin/exp_nagx-cheetah"

echo "=========================================="
echo "exp_nagx Benchmark - VOLE Protocol"
echo "=========================================="
printf "%-17s | %-10s | %-10s | %-19s | %-17s\n" \
    "Elements" "Alice Time" "Bob Time" "Total Comm (bytes)" "Throughput (Exp/s)"
echo "------------------|------------|------------|---------------------|-------------------"

for power in {10..18}
do
    DIM=$((2**power))

    # 清理之前的日志
    rm -f /tmp/alice_${DIM}.log /tmp/bob_${DIM}.log

    # 启动 Alice (服务器)
    timeout 300 $BINARY r=1 d=$DIM nt=1 p=$PORT > /tmp/alice_${DIM}.log 2>&1 &
    ALICE_PID=$!

    # 等待服务器启动
    sleep 2

    # 启动 Bob (客户端)
    timeout 300 $BINARY r=2 d=$DIM nt=1 p=$PORT ip=127.0.0.1 > /tmp/bob_${DIM}.log 2>&1
    BOB_EXIT=$?

    # 等待 Alice 完成
    wait $ALICE_PID 2>/dev/null
    ALICE_EXIT=$?

    # 提取数据（不依赖退出码，直接检查日志）
    ALICE_TIME=$(grep "Exp Time" /tmp/alice_${DIM}.log 2>/dev/null | awk '{print $3}')
    ALICE_COMM=$(grep "Exp Bytes Sent" /tmp/alice_${DIM}.log 2>/dev/null | awk '{print $4}')

    BOB_TIME=$(grep "Exp Time" /tmp/bob_${DIM}.log 2>/dev/null | awk '{print $3}')
    BOB_COMM=$(grep "Exp Bytes Sent" /tmp/bob_${DIM}.log 2>/dev/null | awk '{print $4}')

    # 检查数据是否有效
    if [ -z "$ALICE_TIME" ] || [ -z "$BOB_TIME" ] || [ -z "$ALICE_COMM" ] || [ -z "$BOB_COMM" ]; then
        printf "2^%-2d (%8d)   | FAILED - 检查日志: /tmp/alice_${DIM}.log /tmp/bob_${DIM}.log\n" $power $DIM
        continue
    fi

    # 计算总通信量
    TOTAL_COMM=$((ALICE_COMM + BOB_COMM))

    # 计算吞吐量 (使用较大的时间)
    if [ "$ALICE_TIME" -gt "$BOB_TIME" ] 2>/dev/null; then
        MAX_TIME=$ALICE_TIME
    else
        MAX_TIME=$BOB_TIME
    fi

    THROUGHPUT=$(echo "scale=2; $DIM * 1000 / $MAX_TIME" | bc)

    # 输出结果
    printf "2^%-2d (%8d)   | %7s ms | %7s ms | %19d | %17.0f\n" \
        $power $DIM "$ALICE_TIME" "$BOB_TIME" $TOTAL_COMM $THROUGHPUT

    # 清理进程
    pkill -P $ALICE_PID 2>/dev/null

    # 短暂休息，避免端口冲突
    sleep 1
done

echo "=========================================="
echo "测试完成！"
echo "详细日志: /tmp/alice_*.log 和 /tmp/bob_*.log"
