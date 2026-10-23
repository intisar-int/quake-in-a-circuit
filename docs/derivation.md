# Derivation notes

## 1. Mechanical model

Two lumped masses (floors) with story springs k₁ (ground–floor 1) and k₂ (floor 1–floor 2) and viscous dampers c₁, c₂ in parallel with them. x₁, x₂ are absolute displacements from the fixed ground. The earthquake enters as a sinusoidal force on floor 1.

Floor 1:

```
m₁ x₁'' = F(t) − k₁ x₁ − c₁ x₁' − k₂ (x₁ − x₂) − c₂ (x₁' − x₂')
```

Floor 2:

```
m₂ x₂'' = − k₂ (x₂ − x₁) − c₂ (x₂' − x₁')
```

Matrix form `M x'' + C x' + K x = F`:

```
M = [m₁ 0 ; 0 m₂]
C = [c₁+c₂  −c₂ ; −c₂  c₂]
K = [k₁+k₂  −k₂ ; −k₂  k₂]
F = [F(t) ; 0]
```

## 2. State-space form

With z = [x₁, ẋ₁, x₂, ẋ₂]ᵀ:

```
ż₁ = z₂
ż₂ = −(k₁+k₂)/m₁ z₁ − (c₁+c₂)/m₁ z₂ + k₂/m₁ z₃ + c₂/m₁ z₄ + F(t)/m₁
ż₃ = z₄
ż₄ =  k₂/m₂ z₁ + c₂/m₂ z₂ − k₂/m₂ z₃ − c₂/m₂ z₄
```

C = I₄ (all states observed), D = 0. Inter-story drift is |x₁ − x₂|.

For G = 4:

```
A = [    0      1     0      0
     −1050  −1500    50    500
         0      0     0      1
      1000  10000 −1000 −10000 ]
B = [0 ; 0.00025 ; 0 ; 0]
```

Eigenvalues of A (rad/s): −10552, −946.7, −1.001, −0.1. All are real and negative, so the system is stable and overdamped. The two slow modes (time constants of about 1 s and 10 s) explain why the response is still settling after 5 s.

## 3. Electrical analog (force–voltage / Maxwell)

Mapping: F → V, ẋ → i, x → q, m → L, c → R, k → 1/C.

Series mesh equations, with q₁, q₂ the capacitor charges and i_L1 = q̇₁, i_L2 = q̇₂:

```
V(t) = L₁ i_L1' + R₁ i_L1 + q₁/C₁ + R₂ (i_L1 − i_L2) + (q₁ − q₂)/C₂
0    = L₂ i_L2' + R₂ (i_L2 − i_L1) + (q₂ − q₁)/C₂
```

State vector [q₁, i_L1, q₂, i_L2]ᵀ gives the same A and B as the mechanical model when L = m, R = c and 1/C = k.

| Component | Calculation (G = 4) | Value |
|---|---|---:|
| L₁ | m₁ | 4000 H |
| L₂ | m₂ | 200 H |
| R₁ | c₁ | 4 MΩ |
| R₂ | c₂ | 2 MΩ |
| C₁ | 1/k₁ | 0.25 µF |
| C₂ | 1/k₂ | 5 µF |
| V(t) | F(t) | 5000 sin(8πt) V |

## 4. Sensitivity study

Peak inter-story drift is recomputed with the damping coefficients changed, keeping m and k fixed:

| Case | Change | Peak drift change |
|---|---|---:|
| 1 | c₁ +25 % | −20.42 % |
| 2 | c₂ −25 % | +33.47 % |
| 3 | c₁ and c₂ +25 % | −36.38 % |
| 4 | c₁ and c₂ −25 % | +79.34 % |
