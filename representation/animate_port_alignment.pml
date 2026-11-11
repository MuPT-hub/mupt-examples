import movie_fade

reinitialize
set auto_show_selections, false
load alignment.sdf
select mover, id 19-28
set stick_transparency, 0.5, mover
mset 1-100 x30 100-200
movie_fade stick_transparency, 100, 0.5, 130, 0.0, mover
