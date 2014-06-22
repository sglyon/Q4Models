# API

We have three abstract types that each model should sub-type:

- `AbstractModel`
- `AbstractSolution`
- `AbstractSimulation`

## Fields

Each subtype `AbstractModel` should have fields named `Sol` and `Sim` where the values are types of that `Model`'s subtype. For example, the `IFP` model has the following fields

```
Sol::IFPSol  # where IFPSol <: AbstractSolution
Sim::IFPSim  # where IFPSim <: AbstractSimulation
```

The subtypes of `AbstractSimulation` should all have at least fields for the number of periods in the simulation (named `T::Integer`) and the number of individuals in the panel (named `N::Integer`).

## Methods

Each abstract model should define methods for `solve!(m::CompositeModel, ...; ...)` and `simulate!(m::CompositeModel, ...; ...)`, where `CompositeModel` is the subtype of `AbstractModel`.

These methods can take any number of positional and keyword arguments. For an example see the methods for `IFP` in [ifp.jl]().
