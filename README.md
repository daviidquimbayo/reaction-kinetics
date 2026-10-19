# Reaction Kinetics

Chemical kinetics toolkit in Python.

- Integrated rate laws for zero, first and second order reactions, and their half-lives.
- Arrhenius equation: rate constant vs. temperature and activation energy from two measurements.
- Consecutive reactions `A → B → C` solved numerically with Runge-Kutta 4 and compared with the analytical solution.

## Usage

```bash
pip install -r requirements.txt
python main.py --k1 1.0 --k2 0.5 --a0 1.0 --t 10
pytest
```
