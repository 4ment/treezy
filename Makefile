format:
	black treezy test
	isort treezy test

lint:
	flake8 treezy test

check:
	black --check treezy test
	isort --check-only treezy test
	flake8 treezy test

bump:
	bumpver update
