---
title: ABC444G Solution
date: 2026-02-13 11:17:39
tags:
---

G - Kyoen 题解

> 赛时成功没切 `C` 题，成功被 **8,459** 人吊打

### 分析

题目要求计算圆心为 $(\frac{A}{C}, \frac{B}{C})$、半径为 $\frac{\sqrt{N}}{C}$ 的圆周上格点的数量。

原始方程为：$(x - \frac{A}{C})^2 + (y - \frac{B}{C})^2 = \frac{N}{C^2}$

两边同乘 $C^2$ 得：$(Cx - A)^2 + (Cy - B)^2 = N$

设 $X = Cx - A, Y = Cy - B$，则原问题等价于求解 $X^2 + Y^2 = N$ 满足以下约束的整数解的个数：

- $X \equiv -A \pmod{C}$
- $Y \equiv -B \pmod{C}$

根据费马关于两平方和的定理以及高斯整数环 $\mathbb{Z}[i]$ 的性质：

1. 如果质数 $p \equiv 3 \pmod{4}$ 且在 $N$ 的质因数分解中出现奇数次，则方程 $X^2 + Y^2 = N$ 无整数解。

2. 在 $\mathbb{Z}[i]$ 中，质数的分解规律：
   - $p = 2 = -i(1+i)^2$，其中 $1+i$ 是高斯素数
   - $p \equiv 1 \pmod{4}$ 可分解为 $p = \pi \cdot \bar{\pi}$，其中 $\pi, \bar{\pi}$ 互为共轭高斯素数
   - $p \equiv 3 \pmod{4}$ 在 $\mathbb{Z}[i]$ 中仍为素数

### 思路

首先遍历 $N$ 的质因数分解，如果存在 $p \equiv 3 \pmod{4}$ 且对应指数为奇数，则直接返回 0。

当 $C = 1$ 时，问题退化为标准的两平方和计数问题：答案为 $4 \prod_{p \equiv 1 \pmod{4}} (e_p + 1)$，其中 $e_p$ 是质数 $p$ 在 $N$ 中的指数。

当 $C > 1$ 时，需要考虑同余约束，具体见做法。

### 做法

定义复数运算函数用于处理高斯整数，使用 $cur_{x \cdot C + y}$ 表示当前状态下 $X \equiv x \pmod{C}, Y \equiv y \pmod{C}$ 的方案数。

```cpp
using pii = pair<int, int>;

pii mul(pii a, pii b, int m) {
    ll x = (ll)a.first * b.first - (ll)a.second * b.second;
    x %= m;
    if (x < 0) x += m;
    ll y = (ll)a.first * b.second + (ll)a.second * b.first;
    y %= m;
    if (y < 0) y += m;
    return { (int)x, (int)y };
}

pii pow(pii b, ll e, int m) {
    pii r = {1, 0};
    while (e) {
        if (e & 1) r = mul(r, b, m);
        b = mul(b, b, m);
        e >>= 1;
    }
    return r;
}
```

对于形如 $p \equiv 1 \pmod{4}$ 的质数 $p$，预处理出 $p = a^2 + b^2$ 的表示：

```cpp
vector<int> ps = {5,13,17,29,37,41,53,61,73,89,97};
for (int p : ps) {
    for (int a = 1; a * a <= p; a++) {
        int bb = p - a * a;
        int b = sqrt(bb);
        if (b * b == bb) {
            ra[p] = a;
            rb[p] = b;
            break;
        }
    }
}
```

对于每个质因数 $p$：

- 当 $p \equiv 3 \pmod{4}$：该因子必须成对出现，贡献为 $(p \bmod C)^{e/2}$

```cpp
if (p % 4 == 3) {
    ll f = e / 2;
    int qf = pow_mod(p % C, f, C);
    pii fct = {qf, 0};
    D = mul(D, fct, C);
}
```

- 当 $p \equiv 1 \pmod{4}$：利用 $p = \pi \cdot \bar{\pi}$ 的分解，通过循环卷积计算贡献

```cpp
int a = ra[p];
int b = rb[p];
pii pi = { a % C, b % C };
pii pb = { a % C, ((-b) % C + C) % C };
```

目标状态为 $X \equiv -A \pmod{C}, Y \equiv -B \pmod{C}$，考虑到单位根的影响，需要检查四个可能的目标状态：

- $(X,Y)$, $(−Y,X)$, $(−X,−Y)$, $(Y,−X)$

```cpp
int X0 = (C - A) % C;
int Y0 = (C - B) % C;
int tar = X0 * C + Y0;

vector<int> uid(4);
uid[0] = id(1, 0);
uid[1] = id(0, 1);
uid[2] = id((C - 1) % C, 0);
uid[3] = id(0, (C - 1) % C);

ll ans = 0;
for (int i = 0; i < S; i++) {
    if (cur[i] == 0) continue;
    for (int t = 0; t < 4; t++) {
        int k = mt[i * S + uid[t]];
        if (k == tar) {
            ans = (ans + cur[i]) % MOD;
        }
    }
}
```

最终答案为所有满足目标同余条件的状态的方案数之和，即 $ans$。

复杂度 $O(C^4 + \sum e_i)$，其中 $C$ 为输入参数，$e_i$ 为质因数指数。

由于 DP 部分过于冗长所以解释放代码里了

```cpp
// DP部分：计算满足同余条件的方案数
int S = C * C;  // 状态总数，因为X ≡ x (mod C), Y ≡ y (mod C)，所以共有C*C种状态
int mt[2501][2501]; // 乘法表，mt[i][j]表示状态i和状态j相乘后的结果状态
// 预计算高斯整数乘法的结果，用于快速状态转移
for (int i = 0; i < S; i++) {
    int x1 = i / C;  // 状态i对应的X坐标模C的值
    int y1 = i % C;  // 状态i对应的Y坐标模C的值
    for (int j = 0; j < S; j++) {
        int x2 = j / C;  // 状态j对应的X坐标模C的值
        int y2 = j % C;  // 状态j对应的Y坐标模C的值
        // 计算 (x1 + y1*i) * (x2 + y2*i) = (x1*x2 - y1*y2) + (x1*y2 + x2*y1)*i
        int xr = (x1 * x2 - y1 * y2) % C;  // 实部
        if (xr < 0) xr += C;
        int yr = (x1 * y2 + y1 * x2) % C;  // 虚部
        if (yr < 0) yr += C;
        int k = xr * C + yr;  // 结果状态的编号
        mt[i][j] = k;
    }
}

// 状态编码函数：将坐标(x,y)映射到一维索引
auto id = [&](int x, int y) -> int {
    return x * C + y;
};

// cur[i] 表示当前处理完某些质因数后，处于状态i的方案数
ll cur[2501];
for (int i = 0; i < S; i++) cur[i] = 0;
// 初始化：根据初始因子D设置起始状态
cur[id(D.first, D.second)] = 1;

for (auto& kv : mp) {
    int p = kv.first;   // 当前质数
    ll e = kv.second;   // 当前质数的指数
    if (p % 4 != 1) continue;  // 只处理p ≡ 1 (mod 4)的质数

    // 获取p = a^2 + b^2的分解
    int a = ra[p];
    int b = rb[p];
    // 将质数p在高斯整数环中分解为π = a + bi 和 π̄ = a - bi
    G pi = { a % C, b % C };      // π = a + bi
    G pb = { a % C, ((-b) % C + C) % C };  // π̄ = a - bi

    if (C % p != 0) {  // C与p互质的情况
        // 计算π/π̄的值，用于后续循环节计算
        int nv = (pi.first * pi.first + pi.second * pi.second) % C;  // |π|^2 = π*π̄
        int invn = inv(nv, C);  // |π|^2的逆元
        // 计算π/π̄ = π * π̄ / |π|^2
        G invpb = { (int)((ll)pi.first * invn % C), (int)((ll)pi.second * invn % C) };
        G r = mul(pi, invpb, C);    // r = π/π̄
        G base = pow(pb, e, C);     // base = (π̄)^e

        // 计算循环节长度：找到最小的lam使得r^lam ≡ 1 (mod C)
        int lam = 1;
        G cr = r;
        while (!(cr.first == 1 && cr.second == 0)) {
            cr = mul(cr, r, C);
            lam++;
            if (lam > S) break;  // 防止无限循环
        }

        G seq[2501]; // 存储循环序列中的元素
        G cp = {1, 0};  // 当前幂次
        // 生成序列：base * r^0, base * r^1, ..., base * r^(lam-1)
        for (int t = 0; t < lam; t++) {
            seq[t] = mul(base, cp, C);
            cp = mul(cp, r, C);
        }

        // 统计每个状态在循环序列中出现的次数
        ll cnt[2501];
        for (int i = 0; i < S; i++) cnt[i] = 0;
        for (int t = 0; t < lam; t++) {
            int ix = id(seq[t].first, seq[t].second);
            cnt[ix]++;
        }

        // 计算总项数和循环次数
        unsigned long long tot = e + 1;  // 总共有e+1项需要累加
        unsigned long long full = tot / lam;  // 完整循环的次数
        int rem = tot % lam;  // 剩余不足一次完整循环的项数

        // 构建转移函数f，表示当前质因数对状态转移的贡献
        ll f[2501];
        for (int i = 0; i < S; i++) f[i] = 0;
        ll fm = full % MOD;
        // 对于完整循环的部分
        for (int i = 0; i < S; i++) {
            if (cnt[i] == 0) continue;
            f[i] = fm * (cnt[i] % MOD) % MOD;
        }
        // 对于不完整循环的剩余部分
        for (int t = 0; t < rem; t++) {
            int ix = id(seq[t].first, seq[t].second);
            f[ix] = (f[ix] + 1) % MOD;
        }

        // 执行状态转移：cur = cur * f
        ll nxt[2501];
        for (int i = 0; i < S; i++) nxt[i] = 0;
        for (int i = 0; i < S; i++) {
            if (cur[i] == 0) continue;
            ll ci = cur[i];
            for (int j = 0; j < S; j++) {
                if (f[j] == 0) continue;
                int k = mt[i][j];  // 状态i和状态j相乘后的结果状态
                nxt[k] = (nxt[k] + ci * f[j]) % MOD;
            }
        }
        for (int i = 0; i < S; i++) cur[i] = nxt[i];  // 更新状态
    } else {  // C包含因子p的情况，即p|C
        // 计算C中包含多少个因子p
        int pc = 1;   // p^ce
        int ce = 0;   // p在C中的指数
        int tmp = C;
        while (tmp % p == 0) {
            ce++;
            pc *= p;
            tmp /= p;
        }
        int m = C / pc;  // C = p^ce * m，且gcd(p, m) = 1
        
        ll f[2501];  // 当前质因数的贡献函数
        for (int i = 0; i < S; i++) f[i] = 0;

        if (e < 2LL * ce) {  // 指数较小的情况
            // 直接枚举所有可能的k值：π^k * π̄^(e-k)
            for (ll k = 0; k <= e; k++) {
                G t1 = pow(pi, k, C);      // π^k
                G t2 = pow(pb, e - k, C);  // π̄^(e-k)
                G w = mul(t1, t2, C);      // π^k * π̄^(e-k)
                int ix = id(w.first, w.second);
                f[ix] = (f[ix] + 1) % MOD;
            }
        } else {  // 指数较大的情况，需要优化
            // k < ce 或 k > e-ce 的情况，直接计算
            for (int k = 0; k < ce; k++) {
                G t1 = pow(pi, k, C);
                G t2 = pow(pb, e - k, C);
                G w = mul(t1, t2, C);
                int ix = id(w.first, w.second);
                f[ix] = (f[ix] + 1) % MOD;
            }
            for (ll k = e - ce + 1; k <= e; k++) {
                G t1 = pow(pi, k, C);
                G t2 = pow(pb, e - k, C);
                G w = mul(t1, t2, C);
                int ix = id(w.first, w.second);
                f[ix] = (f[ix] + 1) % MOD;
            }

            // 中间部分：ce <= k <= e-ce，这些项会合并到一起
            ll kl = ce;
            ll kh = e - ce;
            unsigned long long len = (unsigned long long)(kh - kl + 1);

            if (m == 1) {  // 如果C = p^ce，即C只包含p因子
                int zid = id(0, 0);  // 零状态
                f[zid] = (f[zid] + (len % MOD)) % MOD;  // 所有中间项都变为零状态
            } else {  // 如果m > 1，需要更复杂的处理
                // 定义约简函数，将模C降为模m
                auto red = [&](G g) -> G {
                    return { g.first % m, g.second % m };
                };
                G pim = red(pi);  // π模m
                G pbm = red(pb);  // π̄模m

                // 计算π/π̄模m的值
                int nm = (pim.first * pim.first + pim.second * pim.second) % m;
                int invnm = inv(nm, m);
                G invpbm = { (int)((ll)pim.first * invnm % m), (int)((ll)pim.second * invnm % m) };
                G rm = mul(pim, invpbm, m);      // r = π/π̄ (mod m)
                G basem = pow(pbm, e, m);        // base = π̄^e (mod m)

                // 计算模m意义下的循环节长度
                int lam = 1;
                G cr = rm;
                while (!(cr.first == 1 && cr.second == 0)) {
                    cr = mul(cr, rm, m);
                    lam++;
                    if (lam > m * m) break;
                }

                G seqm[2501];  // 模m意义下的循环序列
                G cp = {1, 0};
                for (int t = 0; t < lam; t++) {
                    seqm[t] = mul(basem, cp, m);
                    cp = mul(cp, rm, m);
                }

                int Sm = m * m;
                ll cntm[62501];  // 模m意义下各状态的出现次数
                for (int i = 0; i < Sm; i++) cntm[i] = 0;
                auto idm = [&](int x, int y) -> int {
                    return x * m + y;
                };
                for (int t = 0; t < lam; t++) {
                    int ix = idm(seqm[t].first, seqm[t].second);
                    cntm[ix]++;
                }

                // 计算中间部分的贡献
                unsigned long long full = len / lam;
                int rem = len % lam;
                int off = kl % lam;  // 偏移量

                ll cm[62501];  // 模m意义下的贡献函数
                for (int i = 0; i < Sm; i++) cm[i] = 0;
                ll fm = full % MOD;
                for (int i = 0; i < Sm; i++) {
                    cm[i] = fm * (cntm[i] % MOD) % MOD;
                }
                for (int i = 0; i < rem; i++) {
                    int t = (off + i) % lam;
                    int ix = idm(seqm[t].first, seqm[t].second);
                    cm[ix] = (cm[ix] + 1) % MOD;
                }

                // 将模m的结果提升回模C
                int invpc = inv(pc % m, m);  // (p^ce)^(-1) mod m
                int coord[51];  // 坐标映射数组
                for (int v = 0; v < m; v++) {
                    int tv = (int)((ll)v * invpc % m);  // v * (p^ce)^(-1) mod m
                    int cd = (int)((ll)pc * tv % C);    // 提升回模C
                    coord[v] = cd;
                }

                // 将模m的结果映射到模C
                for (int xm = 0; xm < m; xm++) {
                    for (int ym = 0; ym < m; ym++) {
                        int idxm = xm * m + ym;  // 模m的状态编号
                        ll cntv = cm[idxm];      // 该状态的出现次数
                        if (cntv == 0) continue;
                        int xC = coord[xm];  // 提升实部到模C
                        int yC = coord[ym];  // 提升虚部到模C
                        int idC = xC * C + yC;  // 模C的状态编号
                        f[idC] = (f[idC] + cntv) % MOD;
                    }
                }
            }
        }

        // 执行状态转移：cur = cur * f
        ll nxt[2501];
        for (int i = 0; i < S; i++) nxt[i] = 0;
        for (int i = 0; i < S; i++) {
            if (cur[i] == 0) continue;
            ll ci = cur[i];
            for (int j = 0; j < S; j++) {
                if (f[j] == 0) continue;
                int k = mt[i][j];  // 状态i和状态j相乘后的结果状态
                nxt[k] = (nxt[k] + ci * f[j]) % MOD;
            }
        }
        for (int i = 0; i < S; i++) cur[i] = nxt[i];  // 更新状态
    }
}
```

### [Submission](https://atcoder.jp/contests/abc444/submissions/73128924)
