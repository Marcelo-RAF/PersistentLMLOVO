# PersistentLMLOVO

This project focuses on solving detection problems using an algorithm called **Persistent Levenberg-Marquardt LOVO**.

## Installation Instructions

This package is not yet available in `Metadata.jl`, so it cannot be installed directly using the Julia package manager. To install it, follow these steps:

1. Open the Julia REPL.

2. Press `]` to enter Pkg mode.

3. Run the following command:

   ```julia
   pkg> add https://github.com/Marcelo-RAF/PersistentLMLOVO.git
   ```

4. After installation, press Backspace to return to the Julia REPL.

5. To start using the package, simply run:

   ```julia
   julia> using PersistentLMLOVO
   ```

That's it! Now you can start using the package.

## Running the Algorithm

You can download test examples from the following repository: [PersistentLMLOVO Test Problems](https://github.com/Marcelo-RAF/PersistentLMLOVO/tree/main/problems). Available datasets include:

- **Lines**
- **2D Spheres**
- **3D Spheres**
- **Cubic Functions**
- **Quadratic Functions**

Each category contains problems with only outliers as well as problems with noise and outliers combined.

### Running an Example

To run an example, navigate to the directory where your test data is stored. Suppose you have a file named `line2d_2.0_-4.0_15.csv`, you can load it with:

```julia
julia> prob = load_problem("line2d_2.0_-4.0_15.csv")
```

Once the problem is loaded, you can perform detection using either the **Persistent LMLOVO** or **LMLOVO** algorithm as follows:

```julia
julia> MinimalLMPersistent(xk, prob.model, prob.data, prob.dim, prob.nout)
```

```julia
julia> MinimalLMLOVO(xk, prob.model, prob.data, prob.dim, prob.nout)
```

Where `xk` is the initial guess with dimensions compatible with the problem being analyzed.



### Authors

This project was developed by:

- **Marcelo Renan Augusto Ferreira**
- **Emerson Vitor Castelani**
- **Wesley Vagner Inês Shirabayashi**

