# opengm_rts
opengm_rts is the fast randomness test suits following the GM/T 0005-2021 Randomness test specification.

A lightweight command line executable program is provided to test for 1000 samples with one million bits each.

opengm_rts是随机数检测函数库, 遵循GM/T 0005-2021 随机性检测规范.

opengm_rts可以在10秒内完成对1000组,每组100万比特测试数据的测试.

包括15个随机性测试函数以及命令行程序,方便对1000组1百万比特的样本进行随机性测试, 应用程序也可通过定制api调用实现开机检测和周期检测.

图形界面程序见https://github.com/opengm-libs/opengm_rts_gui

# Build & Usage
Build the command line executable:
```
cargo build --bin opengm_rts --release
```

The command line executable program:
```
$ ./opengm_rts <dir/to/samples>
```

# Performance
1000组,每组100万比特测试数据, Macbook M1 Max上多线程10秒完成.
![performace](/perf.jpg)

# About
OpenGM is the (ongoing) open-source project for the GM/T standards.