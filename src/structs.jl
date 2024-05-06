"""
Structures used throughout package
"""


"""
Mutable structure for options
"""
@with_kw mutable struct EcoNNetOptions
    N::Int64 = 100000;
    window::Int64 = 5000;
    NNet::Bool = true;
    init_weights = glorot_uniform;
    optim = ADAM();
    activation = σ;
    hidden_layers::Int64 = 1;
    train_split::Array{Float64,1} = [1.0, 0.0, 0.0];
    num_nodes::Int64 = 24;
    max_iter::Int64 = 10;
    infoset::Array{Symbol,1};
    expectations::Array{Symbol,1};
    outputs::Array{Symbol,1} = Symbol.([]);
    endogenous::Array{Symbol,1};
    exogenous::Array{Symbol,1};
    states::Array{Symbol,1};
    auxiliary::Array{Symbol,1} = Symbol.([]);
	burnin::Int64 = 6000;
	burnin_use_net::Bool = false;
	learning_gap::Int64 = 100;
	plotting_gap::Int64 = 100;
	plotting_window::Int64 = 999;
	show_plots::Bool = true;
	plot_vars::Array{Symbol,1} = Symbol.([]);
end


"""
Mutable struct with indices of the different types of variable
"""
@with_kw mutable struct EcoNNetIndices
	constant::Int64 = 0;
	infoindex_lagged::Array{Int64,1};
	infonames_lagged::Array{Symbol,1};
	infoindex_current::Array{Int64,1};
	infonames_current::Array{Symbol,1};
	outputindex_current::Array{Int64,1};
	outputnames_current::Array{Symbol,1};
	outputindex_lead::Array{Int64,1};
	outputnames_lead::Array{Symbol,1};
	expectindex_current::Array{Int64,1};
	expectnames_current::Array{Symbol,1};
	expectindex_lead::Array{Int64,1};
	expectnames_lead::Array{Symbol,1};
	expectindex_all::Array{Int64,1};
	expectnames_all::Array{Symbol,1};
	stateindex_lagged::Array{Int64,1};
	statenames_lagged::Array{Symbol,1};
	stateindex_current::Array{Int64,1} ;
	statenames_current::Array{Symbol,1};
	stateindex_all::Array{Int64,1};
	statenames_all::Array{Symbol,1};
	endogindex::Array{Int64,1};
	endognames::Array{Symbol,1};
	auxindex::Array{Int64,1};
	auxnames::Array{Symbol,1};
end
