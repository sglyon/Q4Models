# API

We have three abstract types that each model should sub-type:

- `AbstractModel`
- `AbstractSolution`
- `AbstractSimulation`

Methods defined by API:

- `Base.show`: fully functional show method for printing model to console
- `solve!` and `simulate!`: placeholder methods that raise an error and informs the user that they should define their own version of these methods for their composite subtype. See below for more info.

## Fields

### `AbstractModel`

Each subtype `AbstractModel` **must** have fields named `Sol` and `Sim` where the values are types of that `Model`'s subtype. For example, the `IFP` model has the following fields:

```
...
Sol::IFPSol  # where IFPSol <: AbstractSolution
Sim::IFPSim  # where IFPSim <: AbstractSimulation
...
```

Additionally, if you supply your subtype of `AbstractModel` with a field named `mod_name::String` -- the field value will be used as the title of the model printout.

### `AbstractSimulation`

The subtypes of `AbstractSimulation` **must** all have at least fields for the number of periods in the simulation (named `T::Integer`) and the number of individuals in the panel (named `N::Integer`).

## Methods

Each abstract model should define methods for `solve!(m::CompositeModel, ...; ...)` and `simulate!(m::CompositeModel, ...; ...)`, where `CompositeModel` is the subtype of `AbstractModel`. These methods can take any number of positional and keyword arguments. For an example see the methods for `IFP` in [ifp.jl]().
