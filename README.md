# Symmetric orbits in the restricted three-body problem

This repository contains code for the paper "Computational symplectic topology and symmetric orbits in the restricted three-body problem". The compiled programs give computer assisted proofs of properties of symmetric orbits in the Levi-Civita regularized restricted three-body problem.

## Requirements

The program uses the [CAPD library](http://capd.ii.uj.edu.pl/) for validated numerical computations.

- We used version 5.3.0 of the library, which can be downloaded from [SourceForge](https://sourceforge.net/projects/capd/files/).

- Check out instructions for installation in the [documentation](https://capd.sourceforge.net/capdDynSys/docs/html/a05217.html).

## Compilation
To compile each file execute the following command

```
g++ <file_name>.cpp -o <file_name> `CAPD_BUILD_DIR/bin/capd-config --cflags --libs`
```

where `CAPD_BUILD_DIR` is the directory where CAPD library was installed.

e.g.
```
g++ orbit_properties.cpp -o orbit_properties `~/capd-5.3.0/bin/capd-config --cflags --libs`
g++ retrograde_uniqueness.cpp -o retrograde_uniqueness `~/capd-5.3.0/bin/capd-config --cflags --libs`
g++ symmetric_uniqueness.cpp -o symmetric_uniqueness `~/capd-5.3.0/bin/capd-config --cflags --libs`
g++ retrograde_minimal_action.cpp -o retrograde_minimal_action `~/capd-5.3.0/bin/capd-config --cflags --libs`
g++ direct_minimal_action.cpp -o direct_minimal_action `~/capd-5.3.0/bin/capd-config --cflags --libs`
g++ double_cover_index.cpp -o double_cover_index `~/capd-5.3.0/bin/capd-config --cflags --libs`
```

## Description

The below describes which source code the proof of each theorem relies on.

* **Theorem 1.1**
    - `orbit_properties.cpp`: for existence, non-degeneracy and direct property
    - For example, the following proves the theorem for $\mu \in [0.1, 0.1001]$ and $c \in [1.9, 1.9001]$.
        ```
        ./orbit_properties 1.9 1.9001 0.1 0.1001 -1 0.0 0.001 1 1 64
        ```

* **Theorem 3.1**
    - `orbit_properties.cpp`: for non-degeneracy, retrograde and direct properties
    - `retrograde_uniqueness.cpp`: for uniqueness of retrograde orbit with crossing number 2 
    - `symmetric_uniqueness.cpp`: for uniqueness of symmetric orbits with crossing number 2
    - For example, the following proves the theorem for $\mu \in [0.8, 0.80001]$ and $c \in [2.1, 2.10001]$.
        ```
        ./orbit_properties 2.1 2.10001 0.8 0.80001 1 0.0 0.0001 2 1 8
        ./orbit_properties 2.1 2.10001 0.8 0.80001 -1 0.0 0.0001 2 1 8
        ./symmetric_uniqueness 2.1 2.10001 0.8 0.80001 0 0.0002
        ./symmetric_uniqueness 2.1 2.10001 0.8 0.80001 1 0.0002
        ```
    - As another example, the following proves the theorem for $\mu \in [0.0, 0.0001]$ and $c \in [1.6, 1.6001]$.
        ```
        ./orbit_properties 1.6 1.6001 0.0 0.0001 1 0.0 0.001 2 1 4
        ./orbit_properties 1.6 1.6001 0.0 0.0001 -1 0.5527176322194811 0.01 1 32 0
        ./retrograde_uniqueness 1.6 1.6001 0.0 0.0001 0 0.001
        ./retrograde_uniqueness 1.6 1.6001 0.0 0.0001 1 0.001
        ```

* **Theorem 3.6** (Uniqueness of retrograde orbit as a symmetric orbit)
    - `retrograde_minimal_action.cpp`: for minimal action
    - For example, the following proves the theorem for $\mu \in [0.49998, 0.5]$ and $c \in [2.1, 2.100001]$.
        ```
        ./retrograde_minimal_action 2.1 2.100001 0.49998 0.5 0 0.002
        ./retrograde_minimal_action 2.1 2.100001 0.49998 0.5 1 0.002
        ```

* **Theorem 3.7** (Uniqueness of small action symmetric orbits)
    - `direct_minimal_action.cpp`: for minimal action
    - For example, the following proves the theorem for $\mu \in [0.49998, 0.5]$ and $c \in [2.1, 2.100001]$.
        ```
        ./direct_minimal_action 2.1 2.100001 0.49998 0.5 0 0.0002
        ./direct_minimal_action 2.1 2.100001 0.49998 0.5 1 0.0002
        ```

* **Theorem 3.8**
    - `orbit_properties.cpp`: for non-degeneracy of $\gamma_c$
    - `double_cover_index.cpp`: for index and non-degeneracy of $\gamma_c^2$
    - For example, the following proves non-degeneracy of $\gamma_c$ for $c \in [1.57617, 1.576170005]$.
        ```
        ./orbit_properties 1.57617 1.576170005 0.99 0.99 -1 0.09776738685009109 0.0000005 1 1 64
        ```
    - The index and non-degeneracy of $\gamma_c^2$ are proved by the following
        ```
        ./double_cover_index 1.57617 0.09776738685009109 0.000000000001
        ./double_cover_index 1.57626 0.09965195729775445 0.000000000001
        ```