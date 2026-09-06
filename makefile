all: astral caster sister trial alignment_wtrial chunk_wtrial alignment_wcaster chunk_wcaster

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

trial: dir
	g++ -std=c++20 -march=native -Ofast -D TRIAL src/driver.cpp -o bin/trial

alignment_wtrial: dir
	g++ -std=c++20 -march=native -Ofast -D ALIGNMENT_WTRIAL src/driver.cpp -o bin/alignment_wtrial

chunk_wtrial: dir
	g++ -std=c++20 -march=native -Ofast -D CHUNK_WTRIAL src/driver.cpp -o bin/chunk_wtrial

alignment_wcaster: dir
	g++ -std=c++20 -march=native -Ofast -D ALIGNMENT_WCASTER src/driver.cpp -o bin/alignment_wcaster

chunk_wcaster: dir
	g++ -std=c++20 -march=native -Ofast -D CHUNK_WCASTER src/driver.cpp -o bin/chunk_wcaster

doc: all
	mkdir -p doc
	bin/caster -H > doc/caster.md
