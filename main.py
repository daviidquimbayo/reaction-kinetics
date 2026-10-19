import argparse

import matplotlib.pyplot as plt

from kinetics import consecutive, half_life, intermediate_max_time


def main():
    p = argparse.ArgumentParser(description="Consecutive reactions A → B → C")
    p.add_argument("--k1", type=float, default=1.0, help="rate constant A → B (1/s)")
    p.add_argument("--k2", type=float, default=0.5, help="rate constant B → C (1/s)")
    p.add_argument("--a0", type=float, default=1.0, help="initial [A] (mol/L)")
    p.add_argument("--t", type=float, default=10.0, help="simulation time (s)")
    a = p.parse_args()

    ts, A, B, C = consecutive(a.k1, a.k2, a.a0, t_final=a.t)
    print(f"Half-life of A:     {half_life(1, a.a0, a.k1):.3f} s")
    if a.k1 != a.k2:
        print(f"[B] maximum at t =  {intermediate_max_time(a.k1, a.k2):.3f} s")
    print(f"Final [A], [B], [C]: {A[-1]:.4f}, {B[-1]:.4f}, {C[-1]:.4f} mol/L")

    plt.plot(ts, A, label="[A]")
    plt.plot(ts, B, label="[B]")
    plt.plot(ts, C, label="[C]")
    plt.xlabel("t (s)")
    plt.ylabel("Concentration (mol/L)")
    plt.title(f"A → B → C  (k1 = {a.k1}, k2 = {a.k2})")
    plt.grid(True)
    plt.legend()
    plt.savefig("kinetics.png", dpi=120)
    plt.show()


if __name__ == "__main__":
    main()
