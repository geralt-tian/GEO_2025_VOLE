#!/bin/bash
# SEAF-VOLE-AF exp_nagx测试 - Client端(Bob)
# 使用VOLE (Silent OT) 协议

cd /root/SEAF-VOLE-AF/build/bin

echo "========================================"
echo "启动exp_nagx Client (Bob, party=2)"
echo "协议: VOLE (Silent OT via Ferret)"
echo "========================================"
echo ""

# 参数说明:
# r=2: Bob (client)
# d=1024: 测试1024个元素
# p=32000: 端口号
# ip=127.0.0.1: Server地址（本地测试）
# nt=1: 使用1个线程

./exp_nagx-cheetah r=2 d=1024 p=32000 ip=127.0.0.1 nt=1
