#!/bin/bash
# SEAF-VOLE-AF exp_nagx测试 - Server端(Alice)
# 使用VOLE (Silent OT) 协议

cd /root/SEAF-VOLE-AF/build/bin

echo "========================================"
echo "启动exp_nagx Server (Alice, party=1)"
echo "协议: VOLE (Silent OT via Ferret)"
echo "========================================"
echo ""

# 参数说明:
# r=1: Alice (server)
# d=1024: 测试1024个元素（默认262144太大，先用小数据测试）
# p=32000: 端口号
# nt=1: 使用1个线程

./exp_nagx-cheetah r=1 d=1024 p=32000 nt=1
