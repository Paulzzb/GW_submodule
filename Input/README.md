%%% Purpose
现在的GW module的输入都是直接提供一个struct options_GW，然后修改。
但是实际上这样即使我一个源代码开发者都搞的乱七八糟
现在我期待，重新弄一套规范化输入。
- [ ] 将GW module里面的关键参数分成几大类，
每一大类对应着某一个关键功能。
- [ ] 提供一个稳定的、可扩展的参数读取器。
- [ ] 用markdown书写一个面相用户的参数手册。
- [ ] 测试

Transform the Input of other codes,
such as QE for the ground-state calculation, or BGW for GW calculation,
to the inputs form the code can accept.
