all: astral caster sister caster_tri wtrial

dir:
	g++ -v 2>&1 | tail -n 1
	echo 'If installation failed, please ensure g++ version >= 13'
	mkdir -p bin

astral: dir
	g++ -std=c++20 -march=native -Ofast -D ASTRAL src/driver.cpp -o bin/astral

caster: dir
	g++ -std=c++20 -march=native -Ofast -D CASTER src/driver.cpp -o bin/caster
	
sister: dir
	g++ -std=c++20 -march=native -Ofast -D SISTER src/driver.cpp -o bin/sister

caster_tri: dir
	g++ -std=c++20 -march=native -Ofast -D CASTER_TRI src/driver.cpp -o bin/caster_tri

wtrial: dir
	g++ -std=c++20 -march=native -Ofast -D WEIGHTED_TRIAL src/driver.cpp -o bin/wtrial

doc: all
	mkdir -p doc
	bin/caster -H > doc/caster.md