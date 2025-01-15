- Finish the code in "Input", which transform inputs from other to ours.
- Finish the code that transform the QE result to result we can use.

## 清单
- [ ] 修改Gauss-Lendgre积分公式对应代码.

- [ ] 测试在使用qe的输出结果情形下，计算全频率时候与BGW的差距
  - [ ] 看上去差的不少，一个个修改
- [ ] 现在看上去gpp代码有些问题，需要重新修正一下.
- [ ] 需要 test/ 中添加正确结果，用于进行运行后的对比.
- [ ] 全频率下，用三对角格式生成的 Gauss-Legendre 积分公式, 并测试效果
- [ ] 对多个体系进行奇点测试...
- [ ] A test line
- [ ] 是否需要加一个变量，考虑带k点的兼容（往后放，先只做不带k点的）