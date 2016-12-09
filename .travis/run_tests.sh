#!/bin/bash

coverage run --source=basketball_db --omit="basketball_db/tests" setup.py test
