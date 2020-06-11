<p align="center">
  <img src="https://github.com/mmsbrggr/amber/blob/master/logo.svg">
</p>


# Amber

Amber is an academic prototype to decide the probabilistic termination behavior of Prob-solvalble loops.

## Run Amber with Docker

The repository comes with a Dockerfile. 
This is the easiest way to run Amber.

1. Make sure Docker is installed on your local machine ([docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)).

2. Run the following commands to clone the repository, build the docker image and run the container:

```shell script
git clone git@github.com:mmsbrggr/amber.git
cd amber
docker build -t amber .
docker run -i --name amber1 amber &
```

Now you can run amber on a given benchmark file with the following command:

```shell script
docker exec -t amber1 ./amber benchmarks/past/2d_bounded_random_walk
```

The docker container can be stopped again with executing the following:

```shell script
docker container stop amber1
docker container rm amber1
```


## Local Installation

Amber needs the following dependencies:
- Python version &geq; 3.7 and pip
- scipy
- diofant
- lark-parser

To install these you can do the following steps.

1. Make sure you have python (version &geq; 3.7) and pip installed on your system.
Otherwise install it in your preferred way.

2. Clone the repository:

```shell script
git clone git@github.com:mmsbrggr/amber.git
cd amber
```

3. Create a virtual environment in the `.venv` directory:
```shell script
pip3 install --user virtualenv
python3 -m venv .venv
```

4. Activate the virtual environment:
```shell script
source .venv/bin/activate
```

5. Install the required dependencies with pip:
```shell script
pip install -r requirements.txt
```


## Run Amber

Having all dependencies installed, you can run Amber for example like this:
```shell script
python ./amber.py benchmarks/past/2d_bounded_random_walk
```

A more extensive help can be obtained by:
```shell script
python ./amber.py --help
```

## Writing your own Prob-solvable loop
A Prob-solvable loop consist of initial assignments (one per line), a loop head `while P > Q:`
and a loop body consisting of multiple variable updates (also one per line).
The loop guard `P > Q` consists of two polynomials `P` and `Q` over the program variables.

Initial assignments:
- format:  `var = value
- comment: not all variables have to have initial value specified
- example: `x = 123`

Variable updates:
- format:  `var = option1 @ probability1; option2 @ probability2 ...`
- comment: last probability can be omitted, it's assumed to be whatever
values is needed for probabilities to sum up to 1.
- comment: variables can depend only on previous variables non-linearly,
and on itself linearly - e.g. `x = x + 1` followed by `y = y + x^2` is allowed.
However, `x = x + y` followed by `y = y + 1`, or `x = x^2` is not allowed.
- example: `x = x @ 1/2; x + u`

An example program would be:

```
# this is a comment
y = 0.01
x = 5
while x > 0:
    y = 2*y
    x = x + 200*y**2 @ 1/2; x - y**3
```
More examples can be found in the `benchmarks` folder.
