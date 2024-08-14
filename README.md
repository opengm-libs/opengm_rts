# opengm_rts
opengm_rts is the randomness test suits following the GM/T 0005-2021 Randomness test specification.

A lightweight command line executable program is provided to test for 1000 samples with one million bits each.

opengm_rts是随机数检测函数库, 遵循GM/T 0005-2021 随机性检测规范.

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
For a signle thread test on Apple M1 Max, 1000 samples with one million bits each takes time within 5 minutes:


# OpenGM project
OpenGM is the (ongoing) open-source project for the GM/T standards(SM2, SM3, SM4, SM9, TLCP et al.).