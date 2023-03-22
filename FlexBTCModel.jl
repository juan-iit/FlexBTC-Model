#************************************************************************
# FLEX BTC MODEL
#************************************************************************
# Packages
using JuMP
using Gurobi
using DataFrames
using XLSX

#************************************************************************
# Model
include("Input.jl")
mBTC = Model(Gurobi.Optimizer)

#************************************************************************
# Variables
@variables mBTC begin

    # Electrolizer variables
    m[T]      >= 0      # PV power sold to the grid                    [MWh]
    pv_d[T]   >= 0      # PV power for industrial power demand         [MWh]
    m_d[T]    >= 0      # Power bought from grid for ind power demand  [MWh]
    m_in[T]   >= 0      # Power bought from the grid for standby       [MWh]
    e[T,S]    >= 0      # PV power consumed by electrolyzer by segment [MWh] 
    e_tot[T]  >= 0      # PV power consumed by electrolyzer            [MWh]
    c[T]      >= 0      # PV power consumed by compressor              [MWh]
    h[T]      >= 0      # hydrogen production                          [kg]
    d[T]      >= 0      # hydrogen sold                                [kg]
    s_in[T]   >= 0      # hydrogen stored                              [kg]
    s_out[T]  >= 0      # hydrogen used from storage                   [kg]
    soc[T]    >= 0      # state of charge of storage                   [kg]
    z_h[T,S]  , Bin     # specific hydrogen production                 [-]
    z_on[T]   , Bin     # on electrolyzer                              [-]
    z_off[T]  , Bin     # off electrolyzer                             [-]
    z_sb[T]   , Bin     # standby electrolyzer                         [-]
    z_start[T], Bin     # start electrolyzer                           [-]

    #BTC variables
    p[T]         >= 0   # total BTC power produced                     [MWh]
    p_exp[T]     >= 0   # BTC power sold                               [MWh]
    q[T]         >= 0   # total BTC heat produced                      [MWh]
    q_ns[T]      >= 0   # BTC heat not served                          [MWh]
    q_burn[T]    >= 0   # heat provided by burning biomass             [MWh]
    k_on[T]      , Bin  # on BTC                                       [-]
    k_off[T]     , Bin  # off BTC                                      [-]
    k_start[T]   , Bin  # start BTC                                    [-]

end

#************************************************************************
# Objective function [€]                                                  
# Max profit
@objective(mBTC, Max, 
    sum(  m[t] * lambda_M[t]                  # selling PV power to the grid
        + d[t] * lambda_H                     # selling hydrogen to the market
        + p_exp[t] * lambda_M[t]              # selling BTC power to the grid
        - m_d[t] * lambda_M_in[t]             # cost of purchasing power from the grid to supply the industrial power demand
        - m_in[t] * lambda_M_in[t]            # cost of purchasing power from the grid to supply the electrolyzer during standby
        - z_start[t] * lambda_start           # cost of electrolyzer stratup
        - k_start[t] * D_B * B_price * alpha  # cost of BTC startup if fuelled by biomass
        - q[t] * gamma                        # cost of biomass for BTC heat production
        - p[t] * gamma                        # cost of biomass for BTC power production
        - q_ns[t] * qns_price                 # cost of BTC heat not served
        - q_burn[t] * cf_price                # cost of burning natural gas to supply heat when BTC is OFF
        for t=T))

#************************************************************************
# Constraints

@constraint(mBTC, [t=T],
    m[t]             == P_PV[t] + m_in[t] - e_tot[t] - c[t] - pv_d[t])         # PV power balance

@constraint(mBTC, [t=T],
    pv_d[t] + m_d[t] <= P_demand[t])                                           # Max pv_d and m_d boundary

@constraint(mBTC, [t=T],
    m_in[t]    <= P_sb * z_sb[t])                                              # Limit on standby power from grid

@constraint(mBTC, [t=T],
    e_tot[t]   == sum(e[t,s] for s=S) + P_sb * z_sb[t])                        # Electrolyzer electricity consumption

@constraint(mBTC, [t=T],
    h[t]       == sum(a[s] * e[t,s] + b[s] * z_h[t,s] for s=S))                # Sum segments for hydrogen production 

@constraint(mBTC, [t=T,s=S],
    e[t,s]     >= P_segments[segments][s] * C_E * z_h[t,s])                    # Segment min power boundary

@constraint(mBTC, [t=T,s=1:length(S)-1],
    e[t,s]     <= P_segments[segments][s+1] * C_E * z_h[t,s])                  # Segment max power boundary

@constraint(mBTC, [t=T],
    h[t]       == s_in[t])                                                     # Hydrogen production balance

@constraint(mBTC, [t=T],
    d[t] + (D_H2 * k_start[t] * (1 - alpha)) == s_out[t])                      # Hydrogen output balance

@constraint(mBTC, [t=T],
    s_out[t]   <= C_E * eta_full_load)                                         # Maximum storage output

@constraint(mBTC, [t=T],
    e_tot[t]   <= C_E * z_on[t] + P_sb * z_sb[t])                              # Maximum electrolyzer power

@constraint(mBTC, [t=T],  
    e_tot[t]   >= P_min * z_on[t] + P_sb * z_sb[t])                            # Minimum electrolyzer power

@constraint(mBTC, [t=T],
    z_on[t]    == sum(z_h[t,s] for s=S))                                       # Only one efficiency if on or standby

@constraint(mBTC, [t=T],
    1          == z_on[t] + z_off[t] + z_sb[t])                                # Electrolyzer states

@constraint(mBTC, [t=T],
    1          >= (t > 1 ? z_sb[t] + z_off[t-1] : 0))                          # Not from Off to Standby

@constraint(mBTC, [t=T],
    z_start[t] >= (t > 1 ? z_on[t] - z_on[t-1] - z_sb[t-1] : 0))               # Electrolyzer startup 
 
@constraint(mBTC, [t=T],
    c[t]       == s_in[t] * P_C)                                               # Compressor consumption

@constraint(mBTC, [t=T],
    soc[t]     <= C_S)                                                         # Max storage fill

@constraint(mBTC, [t=T],
    soc[t]     == (t > 1 ? soc[t-1] : 0) - s_out[t] + s_in[t])                 # SOC

@constraint(mBTC, [t=T],
    1          == k_on[t] + k_off[t] + k_start[t])                             # BTC states

# @constraint(mBTC, k_on[1] == 0)

# @constraint(mBTC, k_start[1] == 0)
    
@constraint(mBTC, [t=T],
    1          >= (t > 1 ? k_on[t] + k_off[t-1] : 0))                          # BTC Not from Off to ON
    
@constraint(mBTC, [t=T],
    1          >= (t > 1 ? k_off[t] + k_start[t-1] : 0))                       # BTC Not from Starting to OFF
    
@constraint(mBTC, [t=T],
    1          >= (t > 1 ? k_start[t] + k_on[t-1] : 0))                        # BTC Not from ON to starting

for w in 1:(2 + (alpha*10))
    for t in (w+1):length(T)
        @constraint(mBTC, k_on[t] - k_on[t-1] - k_start[t-w] <= 0)
    end
end

@constraint(mBTC, [t=T],
    p[t]       <= ((m1 * q[t]) + b1) * k_on[t])                                # BTC operation region boundary 1

@constraint(mBTC, [t=T],
    p[t]       >= ((m2 * q[t]) + b2) * k_on[t])                                # BTC operation region boundary 2
 
@constraint(mBTC, [t=T],
    p[t]       >= ((m3 * q[t]) + b3) * k_on[t])                                # BTC operation region boundary 3

@constraint(mBTC, [t=T],
    ((p[t] - p_exp[t]) * k_on[t]) + pv_d[t] + m_d[t] == P_demand[t])           # BTC power balance

@constraint(mBTC, [t=T],
    p_exp[t]   <= (P_no * 1.15) * k_on[t])                                     # Max power_export boundary

@constraint(mBTC, [t=T],
   ((q[t] - q_ns[t]) * k_on[t]) + q_burn[t] == Q_demand[t])                    # BTC heat balance

@constraint(mBTC, [t=T],
    q_ns[t]    <=  Q_no * k_on[t])                                             # Max heat_notserved boundary

@constraint(mBTC, [t=T],
    0          >=  (1 - k_off[t]) * (1 - (P_demand[t] + Q_demand[t])))         # BTC off

#************************************************************************
# Solve
optimize!(mBTC)

# Report results
println("-------------------------------------");
if termination_status(mBTC) == MOI.OPTIMAL
    println(objective_value(mBTC))
    output                    = DataFrame()
    output[!,:H2start]        = value.(z_start) 
    output[!,:H2on]           = value.(z_on)
    output[!,:H2sb]           = value.(z_sb)
    output[!,:H2off]          = value.(z_off)
    output[!,:BTCstart]       = value.(k_start)
    output[!,:BTCon]          = value.(k_on)
    output[!,:BTCoff]         = value.(k_off)
    output[!,:Mprice]         = lambda_M
    output[!,:Mbought]        = value.(m_in)        
    output[!,:PVprod]         = P_PV
    output[!,:PVsold]         = value.(m) 
    output[!,:H2prod_kg]      = value.(h) 
    output[!,:H2stin_kg]      = value.(s_in)
    output[!,:H2stout_kg]     = value.(s_out) 
    output[!,:H2sold_kg]      = value.(d)
    output[!,:SOC_kg]         = value.(soc)
    output[!,:Pdmnd]          = P_demand
    output[!,:Pbtc]           = value.(p)
    output[!,:Psold]          = value.(p_exp)
    output[!,:PVdmnd]         = value.(pv_d)
    output[!,:Mdmd]           = value.(m_d)
    output[!,:Qdmnd]          = Q_demand
    output[!,:Qbtc]           = value.(q)
    output[!,:Qns]            = value.(q_ns)
    output[!,:Qburn]          = value.(q_burn)                         
    show(output) 
    println("\n\nRESULTS:")
    println("Objective value  = $(round(objective_value(mBTC), digits=digs)) EUR")
    println("Electricity sold = $(round(sum((value.(m[t]) + value.(p_exp[t])) * lambda_M[t] for t=T), digits=digs)) EUR")
    println("Hydrogen sold    = $(round(sum(value(d[t]) * lambda_H for t=T), digits=digs)) EUR")
    println("Max SOC          = $(maximum(value.(soc))) kgH2")
    else
    println("No solution")
end
println("\n--------------------------------------");

#************************************************************************
# Output to EXCEL
XLSX.writetable("FlexBTC_Results_C.xlsx", output, overwrite=true)