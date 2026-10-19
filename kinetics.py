"""Chemical kinetics: rate laws, Arrhenius equation and consecutive reactions."""
import math

R = 8.314  # J/(mol·K)


def concentration(order, a0, k, t):
    """[A](t) from the integrated rate law of the given order (0, 1 or 2)."""
    if order == 0:
        return max(a0 - k * t, 0.0)
    if order == 1:
        return a0 * math.exp(-k * t)
    if order == 2:
        return a0 / (1 + a0 * k * t)
    raise ValueError("order must be 0, 1 or 2")


def half_life(order, a0, k):
    if order == 0:
        return a0 / (2 * k)
    if order == 1:
        return math.log(2) / k
    if order == 2:
        return 1 / (k * a0)
    raise ValueError("order must be 0, 1 or 2")


def arrhenius_k(a, ea, temp):
    """k = A·exp(-Ea / RT), with Ea in J/mol and T in K."""
    return a * math.exp(-ea / (R * temp))


def activation_energy(k1, t1, k2, t2):
    """Ea (J/mol) from rate constants measured at two temperatures."""
    return R * math.log(k2 / k1) / (1 / t1 - 1 / t2)


def _rates(a, b, k1, k2):
    return -k1 * a, k1 * a - k2 * b, k2 * b


def consecutive(k1, k2, a0, t_final=10.0, dt=0.001):
    """A → B → C with RK4. Returns (ts, A, B, C)."""
    a, b, c = a0, 0.0, 0.0
    ts, As, Bs, Cs = [0.0], [a], [b], [c]
    n = int(round(t_final / dt))
    for i in range(1, n + 1):
        r1 = _rates(a, b, k1, k2)
        r2 = _rates(a + dt / 2 * r1[0], b + dt / 2 * r1[1], k1, k2)
        r3 = _rates(a + dt / 2 * r2[0], b + dt / 2 * r2[1], k1, k2)
        r4 = _rates(a + dt * r3[0], b + dt * r3[1], k1, k2)
        a += dt / 6 * (r1[0] + 2 * r2[0] + 2 * r3[0] + r4[0])
        b += dt / 6 * (r1[1] + 2 * r2[1] + 2 * r3[1] + r4[1])
        c += dt / 6 * (r1[2] + 2 * r2[2] + 2 * r3[2] + r4[2])
        ts.append(i * dt)
        As.append(a)
        Bs.append(b)
        Cs.append(c)
    return ts, As, Bs, Cs


def intermediate_analytical(k1, k2, a0, t):
    """[B](t) for A → B → C (k1 ≠ k2)."""
    return a0 * k1 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t))


def intermediate_max_time(k1, k2):
    """Time at which [B] reaches its maximum (k1 ≠ k2)."""
    return math.log(k1 / k2) / (k1 - k2)
