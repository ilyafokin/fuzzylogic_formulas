using QuadGK

function gauss(x, A, μ, σ)
    A / sqrt(2 * pi * σ ^ 2) * exp(-0.5 * ((x - μ) / σ) ^ 2)
end

function ∫(f)
    first(quadgk(f, 0., 200.))
end

const N_p_mean = 45
const N_pi_mean = 30
const N_p_sigma = 5
const N_pi_sigma = 3
const ρ_ppi = 0.3

const N_bg_mean = 50
const N_bg_sigma = 10

const σ_p = 5
const σ_pi = 5
const σ_bg = 10

const µ_p = 80
const µ_pi = 75
const μ_bg = 60

const idx_to_str = ["p", "pi", "bg"]

signal = x -> gauss(x, N_p_mean, μ_p, σ_p)
background = x -> gauss(x, N_pi_mean, μ_pi, σ_pi)
bg = x -> gauss(x, N_bg_mean, μ_bg, σ_bg)

function ρ(x, j::String)
    if j === "p"
        signal(x)
    elseif j === "pi"
        background(x)
    elseif j == "bg"
        bg(x)
    else
        throw(ArgumentError("Unknown particle $j"))
    end
end

function ρ(x, j::Integer)
    @assert 0 < j <= length(idx_to_str) "Unknown particle number $j"
    ρ(x, idx_to_str[j])
end

function N(part)
    ∫(x -> ρ(x, part))
end

function ω(x, j)
    # j: particle type
    @assert j in idx_to_str "Unknown particle $j"
    norm = ρ(x, "p") + ρ(x, "pi") + ρ(x, "bg")
    ρ(x, j) / norm
end

function P(x, j)
    norm = ∫(xx -> ρ(xx, j))
    ρ(x, j) / norm
end

function κ1(j, i)
    ∫(x -> ω(x, j) * P(x, i))
end

function κ11(j, k, i)
    term1 = ∫(x -> ω(x, j) * ω(x, k) * P(x, i))
    term2 = ∫(x -> ω(x, j) * P(x, i)) * ∫(x -> ω(x, k) * P(x, i))
    term1 - term2
end

function κ111(j, k, l, i)
    term1 = ∫(x -> ω(x, j) * ω(x, k) * ω(x, l) * P(x, i))
    term2 = ∫(x -> ω(x, j) * ω(x, k) * P(x, i)) * ∫(x -> ω(x, l) * P(x, i))
    term3 = ∫(x -> ω(x, j) * P(x, i)) * ∫(x -> ω(x, k) * ω(x, l) * P(x, i))
    term4 = ∫(x -> ω(x, k) * P(x, i)) * ∫(x -> ω(x, j) * ω(x, l) * P(x, i))
    term5 = ∫(x -> ω(x, j) * P(x, i)) * ∫(x -> ω(x, k) * P(x, i)) * ∫(x -> ω(x, l) * P(x, i))
    term1 - term2 - term3 - term4 + 2 * term5
end

function κ1111(a, b, c, d, i)
    term1 = ∫(x -> ω(x, a) * P(x, i)) * ∫(x -> ω(x, b) * P(x, i)) * ∫(x -> ω(x, c) * P(x, i)) * ∫(x -> ω(x, d) * P(x, i))
    term2 = ∫(x -> ω(x, a) * ω(x, b) * P(x, i)) * ∫(x -> ω(x, c) * P(x, i)) * ∫(x -> ω(x, d) * P(x, i))
    term3 = ∫(x -> ω(x, a) * ω(x, c) * P(x, i)) * ∫(x -> ω(x, b) * P(x, i)) * ∫(x -> ω(x, d) * P(x, i))
    term4 = ∫(x -> ω(x, a) * ω(x, d) * P(x, i)) * ∫(x -> ω(x, b) * P(x, i)) * ∫(x -> ω(x, c) * P(x, i))
    term5 = ∫(x -> ω(x, b) * ω(x, c) * P(x, i)) * ∫(x -> ω(x, a) * P(x, i)) * ∫(x -> ω(x, d) * P(x, i))
    term6 = ∫(x -> ω(x, b) * ω(x, d) * P(x, i)) * ∫(x -> ω(x, a) * P(x, i)) * ∫(x -> ω(x, c) * P(x, i))
    term7 = ∫(x -> ω(x, c) * ω(x, d) * P(x, i)) * ∫(x -> ω(x, a) * P(x, i)) * ∫(x -> ω(x, b) * P(x, i))
    term8 = ∫(x -> ω(x, a) * ω(x, b) * ω(x, c) * P(x, i)) * ∫(x -> ω(x, d) * P(x, i))
    term9 = ∫(x -> ω(x, a) * ω(x, b) * ω(x, d) * P(x, i)) * ∫(x -> ω(x, c) * P(x, i))
    term10 = ∫(x -> ω(x, b) * ω(x, c) * ω(x, d) * P(x, i)) * ∫(x -> ω(x, a) * P(x, i))
    term11 = ∫(x -> ω(x, c) * ω(x, d) * ω(x, a) * P(x, i)) * ∫(x -> ω(x, b) * P(x, i))
    term12 = ∫(x -> ω(x, a) * ω(x, b) * P(x, i)) * ∫(x -> ω(x, c) * ω(x, d) * P(x, i))
    term13 = ∫(x -> ω(x, a) * ω(x, c) * P(x, i)) * ∫(x -> ω(x, b) * ω(x, d) * P(x, i))
    term14 = ∫(x -> ω(x, a) * ω(x, d) * P(x, i)) * ∫(x -> ω(x, b) * ω(x, c) * P(x, i))
    term15 = ∫(x -> ω(x, a) * ω(x, b) * ω(x, c) * ω(x, d) * P(x, i))
    -6 * term1 + 2 * term2 + 2 * term3 + 2 * term4 + 2 * term5 + 2 * term6 + 2 * term7 - term8 - term9 - term10 - term11 - term12 - term13 - term14 + term15
end

function κ1(j::Integer, i::Integer)
    j_str = idx_to_str[j + 1]
    i_str = idx_to_str[i + 1]
    κ1(j_str, i_str)
end


function κ11(j::Integer, k::Integer, i::Integer)
    j_str = idx_to_str[j + 1]
    k_str = idx_to_str[k + 1]
    i_str = idx_to_str[i + 1]
    κ11(j_str, k_str, i_str)
end

function κ111(j::Integer, k::Integer, l::Integer, i::Integer)
    j_str = idx_to_str[j + 1]
    k_str = idx_to_str[k + 1]
    l_str = idx_to_str[l + 1]
    i_str = idx_to_str[i + 1]
    κ111(j_str, k_str, l_str, i_str)
end


function κ1111(a::Integer, b::Integer, c::Integer, d::Integer, i::Integer)
    a_str = idx_to_str[a + 1]
    b_str = idx_to_str[b + 1]
    c_str = idx_to_str[c + 1]
    d_str = idx_to_str[d + 1]
    i_str = idx_to_str[i + 1]
    κ1111(a_str, b_str, c_str, d_str, i_str)
end