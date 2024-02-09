install:
	pip install --upgrade pip &&\
		pip install -r requirements.txt

test:
	python -m pytest -vv -rP --cov=hello test_units.py

lint:
	pylint --disable=R,C domain_mapper.py

format:
	black *.py

all: install test
