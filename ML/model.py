# model.py
import torch
import torch.nn as nn
import torch.nn.functional as F

class MLP(nn.Module):
    def __init__(self, in_dim, hidden=[128, 128], act=nn.SiLU):
        super().__init__()
        layers = []
        dims = [in_dim] + hidden
        for i in range(len(dims)-1):
            layers += [nn.Linear(dims[i], dims[i+1]), act()]
        self.body = nn.Sequential(*layers)

    def forward(self, x):
        return self.body(x)

class StabilityGate(nn.Module):
    def __init__(self, init_alpha=5.0, temperature=1.0):
        super().__init__()
        self.log_alpha = nn.Parameter(torch.tensor(init_alpha).log())
        self.log_temp  = nn.Parameter(torch.tensor(temperature).log())

    def forward(self, zeta):
        alpha = self.log_alpha.exp()
        temp  = self.log_temp.exp().clamp(min=1e-2, max=10.0)
        return torch.sigmoid(alpha * zeta / temp)

class FluxNN(nn.Module):
    def __init__(self, in_dim, hidden=[128,128]):
        super().__init__()
        self.gate = StabilityGate()
        self.stable = MLP(in_dim, hidden=hidden)
        self.unstable = MLP(in_dim, hidden=hidden)
        # heads for tau and H
        self.head_tau_stable   = nn.Linear(hidden[-1], 1)
        self.head_tau_unstable = nn.Linear(hidden[-1], 1)
        self.head_H_stable     = nn.Linear(hidden[-1], 1)
        self.head_H_unstable   = nn.Linear(hidden[-1], 1)

    def forward(self, x, zeta):
        # x: features; zeta: requires_grad tensor
        g = self.gate(zeta)
        f_stable   = self.stable(x)
        f_unstable = self.unstable(x)
        tau_st = self.head_tau_stable(f_stable)
        tau_un = self.head_tau_unstable(f_unstable)
        H_st   = self.head_H_stable(f_stable)
        H_un   = self.head_H_unstable(f_unstable)

        tau_hat = g * tau_un + (1 - g) * tau_st
        H_hat   = g * H_un   + (1 - g) * H_st
        return tau_hat.squeeze(-1), H_hat.squeeze(-1), g
