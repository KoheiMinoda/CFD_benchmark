# 角柱まわり流れの数値解析設定

正方形断面角柱（一辺 $H$）まわりの渦放出流れ、$Re = 21{,}400$。
ERCOFTAC Classic Collection case043（Lyn, Einav, Rodi & Park, *J. Fluid Mech.* **304** (1995) 285–319）の
実験条件に合わせた二次元非定常 RANS 解析。

---

## 1. 記号と基準量

| 記号 | 意味 | 値 |
|---|---|---|
| $H$ | 角柱の一辺（代表長さ） | $0.04\ \mathrm{m}$ |
| $U_0$ | 一様流速（代表速度） | $0.535\ \mathrm{m\,s^{-1}}$ |
| $\nu$ | 動粘性係数 | $1.0\times10^{-6}\ \mathrm{m^2\,s^{-1}}$ |
| $\rho$ | 密度 | $1.0\times10^{3}\ \mathrm{kg\,m^{-3}}$ |
| $Re$ | $U_0 H/\nu$ | $21{,}400$ |
| $t^{*}$ | 無次元時間 $tU_0/H$ | — |

座標原点は角柱断面の中心。$x$ は主流方向、$y$ は横断方向、$z$ はスパン方向。

---

## 2. 支配方程式

非圧縮性・一定物性の非定常レイノルズ平均 Navier–Stokes 方程式（URANS）を解く。
$p$ は密度で割った**修正圧力**（後述）である。

$$
\frac{\partial U_i}{\partial x_i} = 0
$$

$$
\frac{\partial U_i}{\partial t}
+ \frac{\partial}{\partial x_j}\left( U_i U_j \right)
= -\frac{\partial p}{\partial x_i}
+ \frac{\partial}{\partial x_j}\left[
\left( \nu + \nu_t \right)
\left( \frac{\partial U_i}{\partial x_j} + \frac{\partial U_j}{\partial x_i} \right)
\right]
$$

レイノルズ応力は Boussinesq 渦粘性近似で閉じる。

$$
-\overline{u_i' u_j'} = \nu_t \left( \frac{\partial U_i}{\partial x_j} + \frac{\partial U_j}{\partial x_i} \right) - \frac{2}{3} k\, \delta_{ij}
$$

**注意（他コードとの比較時に重要）**：等方成分 $\tfrac{2}{3}k\,\delta_{ij}$ は独立に扱わず圧力に繰り込んでいる。
したがって解かれている変数は

$$
p = \frac{P}{\rho} + \frac{2}{3}k
\qquad [\mathrm{m^2\,s^{-2}}]
$$

であり、圧力係数を厳密に比較する場合はこの $\tfrac{2}{3}k$ の差（本ケースでは $2k/3U_0^2 \lesssim 0.05$）を考慮する。

---

## 3. 乱流モデル：$k$–$\omega$ SST（Menter 2003）

### 3.1 輸送方程式

$$
\frac{\partial k}{\partial t} + \frac{\partial (U_j k)}{\partial x_j}
= \tilde{P}_k - \beta^{*} k \omega
+ \frac{\partial}{\partial x_j}\left[ \left( \nu + \sigma_k \nu_t \right) \frac{\partial k}{\partial x_j} \right]
$$

$$
\frac{\partial \omega}{\partial t} + \frac{\partial (U_j \omega)}{\partial x_j}
= \gamma S^2 - \beta \omega^2
+ \frac{\partial}{\partial x_j}\left[ \left( \nu + \sigma_\omega \nu_t \right) \frac{\partial \omega}{\partial x_j} \right]
+ 2 (1 - F_1) \frac{\sigma_{\omega 2}}{\omega} \frac{\partial k}{\partial x_j}\frac{\partial \omega}{\partial x_j}
$$

ここで
$S_{ij} = \tfrac{1}{2}\left( \partial U_i/\partial x_j + \partial U_j/\partial x_i \right)$、
$S^2 = 2 S_{ij}S_{ij}$。

### 3.2 渦粘性と生成項制限

$$
\nu_t = \frac{a_1 k}{\max\!\left( a_1 \omega,\; b_1 F_2 \sqrt{S^2} \right)}
$$

$$
P_k = \nu_t S^2, \qquad
\tilde{P}_k = \min\!\left( P_k,\; c_1 \beta^{*} k \omega \right)
$$

$\omega$ 方程式の生成項にも同じ比率の制限を課している。

### 3.3 混合関数

$$
F_1 = \tanh\!\left( \Gamma_1^4 \right), \qquad
\Gamma_1 = \min\left[ \max\left( \frac{\sqrt{k}}{\beta^{*}\omega d},\; \frac{500\nu}{d^2 \omega} \right),\;
\frac{4 \sigma_{\omega 2} k}{CD_{k\omega} d^2} \right]
$$

$$
CD_{k\omega} = \max\!\left( \frac{2\sigma_{\omega 2}}{\omega} \frac{\partial k}{\partial x_j}\frac{\partial \omega}{\partial x_j},\; 10^{-10} \right)
$$

$$
F_2 = \tanh\!\left( \Gamma_2^2 \right), \qquad
\Gamma_2 = \max\left( \frac{2\sqrt{k}}{\beta^{*}\omega d},\; \frac{500\nu}{d^2 \omega} \right)
$$

$d$ は最近傍壁面までの距離。前進波面法（wave-front / fast-marching 型アルゴリズム）で厳密な最短距離を評価している。

モデル定数は $\phi = F_1 \phi_1 + (1-F_1)\phi_2$ で内層・外層の値をブレンドする。

### 3.4 モデル定数

| 定数 | 値 | 定数 | 値 |
|---|---|---|---|
| $\sigma_{k1}$ | 0.85 | $\sigma_{k2}$ | 1.0 |
| $\sigma_{\omega 1}$ | 0.5 | $\sigma_{\omega 2}$ | 0.856 |
| $\beta_1$ | 0.075 | $\beta_2$ | 0.0828 |
| $\gamma_1$ | 5/9 | $\gamma_2$ | 0.44 |
| $\beta^{*}$ | 0.09 | $\kappa$ | 0.41 |
| $a_1$ | 0.31 | $b_1$ | 1.0 |
| $c_1$ | 10 | | |

いずれも Menter (2003) の標準値。曲率補正・回転補正・遷移モデルは使用していない。

---

## 4. 計算領域

二次元平面（$x$–$y$）。奥行き方向は 1 セル層のみで、両側面に二次元条件を課している。

| 項目 | 値 | 無次元 |
|---|---|---|
| 上流長さ（流入面〜角柱中心） | $0.8\ \mathrm{m}$ | $20H$ |
| 下流長さ（角柱中心〜流出面） | $2.0\ \mathrm{m}$ | $50H$ |
| 流路幅 | $0.56\ \mathrm{m}$ | $14H$（$-7H \le y/H \le 7H$）|
| 奥行き | $0.004\ \mathrm{m}$ | $0.1H$（1 セル）|
| 閉塞率 $H/W$ | — | $7.14\%$ |

流路幅は実験の水路断面（$0.39\times0.56\ \mathrm{m}$ のうち $0.56\ \mathrm{m}$ 側）に一致させてあり、
閉塞率も実験の $7\%$ と揃っている。上下面は実験と同様に**滑りなし壁**である（滑り条件でも対称条件でもない）。

上流 $20H$ は、角柱前方の圧力上昇が流入面の一様条件に届かない距離であり、
$C_p$ の基準断面（$x/H = -10$、§12）を $y$ 方向に一様な位置に取るためにも必要である。
下流 $50H$ は、流出境界の影響を後流統計の評価域（$x/H \le 8$）から遠ざけるために取っている。

---

## 5. 境界条件

| 境界 | $\boldsymbol{U}$ | $p$ | $k$ | $\omega$ |
|---|---|---|---|---|
| 流入 | $\boldsymbol{U} = (U_0, 0, 0)$ 一様 | $\partial p/\partial n = 0$ | 固定値（下式） | 固定値（下式） |
| 流出 | $\partial \boldsymbol{U}/\partial n = 0$ | $p = 0$ | $\partial k/\partial n = 0$ | $\partial \omega/\partial n = 0$ |
| 角柱表面 | $\boldsymbol{U} = 0$ | $\partial p/\partial n = 0$ | 近傍プロファイル則（§6.1） | ブレンド則（§6） |
| 流路上下壁 | $\boldsymbol{U} = 0$ | $\partial p/\partial n = 0$ | $\partial k/\partial n = 0$（§6.2） | ブレンド則（§6） |
| 側面（$z$） | 二次元条件（$\partial/\partial z = 0$, $w = 0$） | 同左 | 同左 | 同左 |

### 流入乱流量

乱れ強度 $I = 2\%$（実験の自由流乱れ強度に一致）と乱流混合長 $\ell = 0.01\ \mathrm{m} = 0.25H$ から

$$
k_\infty = \frac{3}{2}\left( I\, U_0 \right)^2 = 1.717\times10^{-4}\ \mathrm{m^2\,s^{-2}}
\qquad \left( k_\infty / U_0^2 = 6.0\times10^{-4} \right)
$$

$$
\omega_\infty = \frac{\sqrt{k_\infty}}{C_\mu^{1/4}\,\ell} = 2.393\ \mathrm{s^{-1}}
\qquad \left( \omega_\infty H/U_0 = 0.179 \right), \quad C_\mu = \beta^{*} = 0.09
$$

$$
\left. \frac{\nu_t}{\nu} \right|_\infty = \frac{k_\infty}{\omega_\infty \nu} = 71.8
$$

初期条件は全領域で $\boldsymbol{U} = (U_0,0,0)$、$p=0$、$k = 4.8\times10^{-4}\ \mathrm{m^2 s^{-2}}$、$\omega = 1.1\ \mathrm{s^{-1}}$。
統計は初期過渡を十分に過ぎた後で取るため、初期条件は結果に影響しない。

---

## 6. 壁面処理

壁面処理は角柱表面と流路上下壁で異なるが、**渦粘性の壁面値はどちらも $\nu_t|_w = 0$** とする。
すなわち壁せん断応力を対数則から求めるのではなく、粘性応力として直接評価する。

$$
\nu_t \big|_{w} = 0
\quad \Longrightarrow \quad
\tau_w = \rho\,\nu \left. \frac{\partial U_t}{\partial n} \right|_{w}
$$

$\omega$ の壁隣接値も両者共通で、粘性底層の解析解と対数層の解を二乗平均で接続した式で規定する。

$$
\omega_{\mathrm{vis}} = \frac{6\nu}{\beta_1 d^2}, \qquad
\omega_{\mathrm{log}} = \frac{\sqrt{k}}{C_\mu^{1/4} \kappa d}, \qquad
\omega_w = \sqrt{\omega_{\mathrm{vis}}^2 + \omega_{\mathrm{log}}^2}
$$

角柱表面（$y^{+} \lesssim 11$）では $\omega_{\mathrm{vis}}$ が、上下壁（$y^{+} \approx 30$）では
$\omega_{\mathrm{log}}$ が支配的になる。

### 6.1 角柱表面：壁まで解像する低レイノルズ数型

$k$ の壁隣接値は普遍的な近傍プロファイル $k^{+}(y^{+}) = k/u_\tau^2$ から与える。

$$
k^{+} =
\begin{cases}
\dfrac{2400}{C_{\varepsilon 2}^{2}}
\left[ \dfrac{1}{(y^{+}+C)^2} + \dfrac{2y^{+}}{C^3} - \dfrac{1}{C^2} \right]
& y^{+} \le y^{+}_{\mathrm{lam}} \\[2ex]
\dfrac{C_k}{\kappa}\ln y^{+} + B_k & y^{+} > y^{+}_{\mathrm{lam}}
\end{cases}
$$

$C_{\varepsilon 2} = 1.9$, $C_k = -0.416$, $B_k = 8.366$, $C = 11.0$,
$y^{+}_{\mathrm{lam}} \simeq 11.53$（線形則と対数則の交点、$\kappa = 0.41$, $E = 9.8$）。

### 6.2 流路上下壁：$k$ は勾配ゼロ

上下壁は $y^{+} \approx 30$ にあり、粘性底層内を前提とした $k^{+}$ プロファイル則を課すのは
適切でない。ここでは $k$ を壁面で勾配ゼロとする。

$$
\left. \frac{\partial k}{\partial n} \right|_{w} = 0
$$

## 7. 計算格子

境界適合の**構造格子**（全六面体）。角柱を囲む O 型ブロックを外側の直交格子に接続した構成。
O 型ブロックは角柱表面から法線方向に一定距離 $1.5H$ だけ外側に置いた正方形（$|x|,|y| \le 2H$）まで
広がり、その内側で壁面から後流形成域にかけて格子を連続的に粗くしていく。
O 型ブロックの外側は等間隔の直交格子である。

| 項目 | 値 | 無次元 |
|---|---|---|
| 総セル数 | 103,920 | — |
| 節点数 | 209,680 | — |
| 角柱表面のセル数 | 160（1 辺あたり 40） | 接線方向間隔 $H/40 = 0.025H$ |
| 壁隣接セル高さ $\Delta y_1$ | $1.98\times10^{-4}\ \mathrm{m}$ | $0.0050H$（$= H/200$）|
| 角柱法線方向のセル数 | 47（$d \le 1.5H$ の範囲）| セルあたり伸長率 1.067（最大／最小 = 20）|
| 後流・遠方の等間隔セル | $4.0\times10^{-3}\ \mathrm{m}$ | $0.1H$（$= H/10$）|

角柱の四隅では O 型ブロックの稜線が対角方向に走るため、第 1 セルの高さは辺の中央の
$\sqrt{2}$ 倍になる。

### O 型ブロック内の格子幅

| 壁からの距離 $d$ | そこまでのセル数 | 局所格子幅 |
|---|---|---|
| $0.05H$ | 7 | $0.0078H$ |
| $0.10H$ | 13 | $0.0116H$ |
| $0.20H$ | 20 | $0.0182H$ |
| $0.50H$ | 31 | $0.0374H$ |
| $1.00H$ | 41 | $0.0716H$ |
| $1.50H$ | 47 | $0.0992H$ |

後流中心線でいえば、角柱背面から $x/H = 1$ で $\Delta \approx 0.037H$、$x/H = 1.5$ で $0.072H$、
$x/H = 2$ で $0.099H$、それ以遠は $0.1H$ 一定である。せん断層厚さ $0.2$–$0.3H$ に対して
渦形成域（$0.5 \le x/H \le 2$）で 2–8 セル、前縁剥離直後（$d \le 0.1H$）で 13 セルを充てている。

### 格子品質

| 指標 | 値 |
|---|---|
| 最大非直交角 | $44.3^{\circ}$（平均 $7.4^{\circ}$）|
| 最大歪度（skewness） | 2.50 |
| 最大アスペクト比 | 6.04 |
| 最小セル体積 | $7.98\times10^{-10}\ \mathrm{m^3}$ |

---

## 8. 空間離散化

セル中心配置・非スタガード（コロケート）配置の**有限体積法**。すべての項を検査体積で積分し、
面積分は各面の中心値に面積を乗じた二次精度中点則で評価する。

| 項 | 手法 | 精度 |
|---|---|---|
| 勾配 | Green–Gauss セル基準（面値は線形内挿） | 2 次 |
| 運動量の対流項 | 流束制限付き中心差分（下式） | 2 次（極値近傍で 1 次へ低下）|
| $k$, $\omega$ の対流項 | 同上（スカラ版） | 同上 |
| 拡散項 | 中心差分 + 非直交性の陽的補正 | 2 次 |
| 面法線勾配 | 非直交補正付き（over-relaxed approach） | 2 次 |
| その他の面内挿 | 線形（中心） | 2 次 |

### 対流項の流束制限

面値は 1 次風上値と 2 次中心値の混合で与える。

$$
\phi_f = \phi_f^{\mathrm{UD}} + \psi(r)\left( \phi_f^{\mathrm{CD}} - \phi_f^{\mathrm{UD}} \right),
\qquad
\psi(r) = \max\!\left[ \min\left( 2r,\, 1 \right),\, 0 \right]
$$

$r$ は風上側・下流側の勾配比。$\psi = 1$ で純中心差分、$\psi = 0$ で純風上差分。
この制限関数は Sweby 線図上で TVD 条件を満たす（$\beta = 2$ の Osher 型に相当）。

運動量方程式では制限関数を速度ベクトルに対して**一括で**評価する。すなわち変化が最も急な方向で
$r$ を求め、得られた $\psi$ を 3 成分すべてに適用する（成分ごとに別々の制限をかけない）。
これは成分間の不整合による非物理的な変形を避けるための処置である。

### 非直交性補正

面法線勾配は直交成分と非直交補正成分に分解し、後者を陽的に扱う。

$$
\left. \frac{\partial \phi}{\partial n} \right|_f
= \underbrace{\frac{\phi_N - \phi_P}{|\boldsymbol{d}|}\,\frac{|\boldsymbol{S}_f|}{\boldsymbol{S}_f \cdot \hat{\boldsymbol{d}}}}_{\text{陰的}}
+ \underbrace{\left( \boldsymbol{S}_f - \frac{|\boldsymbol{S}_f|^2}{\boldsymbol{S}_f \cdot \boldsymbol{d}} \boldsymbol{d} \right) \cdot \overline{\nabla \phi}_f}_{\text{陽的補正}}
$$

補正は圧力方程式の中で 1 回追加反復して収束させる（§10）。

### 乱流量の対流スキームについて

$k$, $\omega$ の対流項にも運動量と同じ流束制限付き中心差分を用いている。この組み合わせでは
$\omega$ が局所的に負値に振れて下限クリップされる時間ステップが全体の 4.2%、$k$ で 34.7% 生じる
（クリップ量は $\omega$ の代表値に対して小さいが、非保存的な生成項として働く）。
クリッピングを完全に消したい場合は乱流量の対流のみ 1 次風上（`Gauss upwind`）にすればよい。
運動量の対流スキームは 2 次のままなので全体の精度は落ちない。

---

## 9. 時間離散化

**1 次精度陰的 Euler 法**（後退差分）。

$$
\frac{\partial \phi}{\partial t} \approx \frac{\phi^{n+1} - \phi^{n}}{\Delta t}
$$

時間刻みは局所クーラン数の上限で自動調整する。

$$
\mathrm{CFL} = \Delta t \max_{\text{cells}} \frac{\sum_f |\boldsymbol{U}_f \cdot \boldsymbol{S}_f|}{2 V_P} \le 0.8,
\qquad \Delta t \le 1.0\times10^{-3}\ \mathrm{s}
$$

| 項目 | 値 |
|---|---|
| 最大クーラン数 | 0.8 |
| 実効時間刻み | $\Delta t \approx 2.48\times10^{-4}\ \mathrm{s}$（$\Delta t\, U_0/H = 3.31\times10^{-3}$）|
| 総時間ステップ数 | 242,260 |
| 全計算時間 | $60\ \mathrm{s}$（$t^{*} = 802$、渦放出 約 103 周期）|
| 1 放出周期あたりのステップ数 | 約 2,340 |

---

## 10. 圧力–速度連成

分離解法（segregated）による**射影法**。1 時間ステップあたりの手順は次のとおり。

1. 運動量方程式を陰的に解き、予測速度 $\boldsymbol{U}^{*}$ を得る
2. コロケート格子の圧力–速度デカップリングを防ぐため、面流束は運動量方程式に基づく内挿
   （Rhie–Chow 型の運動量内挿）で構成する
3. 圧力 Poisson 方程式を解く

$$
\nabla \cdot \left( \frac{1}{a_P} \nabla p \right)
= \nabla \cdot \left( \frac{\boldsymbol{H}(\boldsymbol{U})}{a_P} \right)
$$

4. 速度と面流束を補正して連続の式を満たす
5. 乱流量 $k$, $\omega$ を解く

| 反復レベル | 回数 |
|---|---|
| 外部反復（非線形項・連成の更新） | 2 |
| 圧力補正の内部反復 | 各外部反復あたり 2 |
| 非直交性の追加補正 | 各圧力解あたり 1 |
| 陽的緩和係数 | 1.0（緩和なし） |

外部反復を 2 回行うため、非定常性は時間刻みではなく反復で担保される
（緩和を用いない純粋な非定常解法）。

---

## 11. 線形方程式ソルバ

| 方程式 | 解法 | 収束判定 |
|---|---|---|
| 圧力 | 代数マルチグリッド（V サイクル、平滑化は不完全 Cholesky 前処理付き Gauss–Seidel）| 絶対残差 $10^{-8}$、相対残差 $10^{-2}$（各時間ステップの最終反復では相対残差の判定を外し絶対残差 $10^{-8}$ まで）|
| 運動量・$k$・$\omega$ | 対称 Gauss–Seidel 反復 | 同上 |

残差はシステム行列の対角要素と場のスケールで正規化した $L_1$ ノルム。

---

## 12. 統計処理と評価量の定義

### 平均操作

初期過渡を除いた $t = 20$–$60\ \mathrm{s}$（渦放出 69 周期）で時間平均を取る。
後流プロファイルは $t = 40$–$60\ \mathrm{s}$（34.5 周期）を $\Delta t^{*}_{\mathrm{sample}} = 0.379$ 間隔で
707 スナップショット抽出して平均した（位相平均は 20 ビンなので 1 ビンあたり約 35 サンプル、
上下対称性で畳んで約 70）。

### 力係数

$$
C_D = \frac{F_x}{\tfrac{1}{2}\rho U_0^2 H L_z},
\qquad
C_L = \frac{F_y}{\tfrac{1}{2}\rho U_0^2 H L_z},
\qquad
L_z = 0.004\ \mathrm{m}
$$

$F_x$, $F_y$ は角柱表面上の圧力力と粘性力の合計。基準面積は $H \times L_z$（奥行き 1 セル分）。

### ストローハル数

$$
St = \frac{f_s H}{U_0}
$$

$f_s$ は $C_L$ 時系列のパワースペクトル密度（Hann 窓、Welch 法、重複率 50%）の主ピーク周波数。
上向きゼロ交差の平均間隔からも独立に算出して整合を確認する。

### 圧力係数

$$
C_p = \frac{p - p_{\mathrm{ref}}}{\tfrac{1}{2} U_0^2}
$$

$p$ は §2 の修正運動学的圧力。基準 $p_{\mathrm{ref}}$ は角柱中心から $10H$ 上流の鉛直断面
$x/H = -10$（$|y/H| \le 2.5$）における時間・空間平均値を用いる。流出面（$p=0$）を基準にすると
上流静圧のオフセットがそのまま乗り、よどみ点で $C_p > 1$ という非物理的な値になる。
本ケースの実測値は $p_{\mathrm{ref}}/U_0^2 = +0.0976$ で、このときよどみ点 $C_p = +0.991$ となり
理論上限 $+1.0$ と整合する（$x/H = -3$ ではまだ角柱の閉塞による圧力上昇が $y$ 方向に残るため
基準断面には使えない）。

### 再循環長さ

後流中心線 $y = 0$ 上で $\overline{U}$ が負から正に転じる位置 $x_R$（線形内挿）。
角柱中心を原点として無次元化する（角柱背面からの長さは $x_R/H - 0.5$）。

### 乱れ量の分解

URANS では速度変動が「解像された周期成分」と「モデル化された成分」に分かれる。
実験の $u'$, $v'$ と比較する際は両者を合成する。

$$
u'_{\mathrm{tot}} = \sqrt{ \overline{\left(U - \overline{U}\right)^2} + \tfrac{2}{3}k },
\qquad
\overline{u'v'}_{\mathrm{tot}} = \overline{\left(U - \overline{U}\right)\left(V - \overline{V}\right)} - \nu_t \frac{\partial \overline{U}}{\partial y}
$$

（せん断応力のモデル成分は鉛直線上で評価するため $\partial V/\partial x$ を無視した近似）

### 位相平均

$C_L$ を放出周波数まわりで帯域通過（4 次 Butterworth、$0.6f_s$–$1.6f_s$、零位相）した後、
Hilbert 変換で瞬時位相 $\varphi(t)$ を求め、$18^{\circ}$ 刻みの 20 ビンにスナップショットを分類する。

$$
\langle \phi \rangle (\boldsymbol{x}, \varphi_m) = \frac{1}{N_m} \sum_{t_i \in \varphi_m} \phi(\boldsymbol{x}, t_i),
\qquad
\langle \phi' \rangle = \sqrt{ \frac{1}{N_m} \sum_{t_i \in \varphi_m} \left( \phi - \langle \phi \rangle \right)^2 }
$$

流れ場の反対称性 $\langle \phi \rangle (x, -y, \varphi) = \pm \langle \phi \rangle (x, y, \varphi + 180^{\circ})$
（$U$, $u'$, $v'$ は偶、$V$, $\overline{u'v'}$ は奇）を用いて上下半面を重ね、標本数を倍にしている。

---

## 参考文献

1. Menter, F. R., Kuntz, M. & Langtry, R. (2003) *Ten years of industrial experience with the SST turbulence model.* Turbulence, Heat and Mass Transfer 4, 625–632.
2. Lyn, D. A., Einav, S., Rodi, W. & Park, J.-H. (1995) *A laser-Doppler velocimetry study of ensemble-averaged characteristics of the turbulent near wake of a square cylinder.* J. Fluid Mech. **304**, 285–319.
3. ERCOFTAC Classic Collection, Case 043: Vortex Shedding Past Square Cylinder. <http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case043>
