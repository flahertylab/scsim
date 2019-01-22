# SCSIM: Jointly Simulating Correlated single-cell and bulk next-generation sequencing data

This is a tool for simulating  next-generation sequencing data from a hierarchical sampling arrangement with single-cell and bulk samples.

## Getting Started

1. Clone the repository to your computer.
2. Inside the repository, run `make build` to build the =scsim= docker image.

### Prerequisites

- docker

### Installing 

## Running the tests

The =scsim= docker image must be built before running the example. 

1. Navigate the the =examples= directory in the repository.
2. Run `make all`

A docker container will be created and the example simulation will be run in the container.
The simulation output will appear in the =example= folder because a shared volume is created in the container that connects to the =example= folder.

## Contributing

## Versioning

## Authors

* Collin Giguere
* Patrick Flaherty

## License

## Acknowledgements
