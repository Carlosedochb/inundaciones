model WaterOnFields

global {
	grid_file dem_file <- file("../includes/gdl_pdz.tif");
	field terrain <- field(dem_file);
	field flow <- field(dem_file);
	geometry shape <- envelope(dem_file);
	
	bool fill <- false;
	bool use_dem_elevation <- true; 
	bool rainfall_mode <- true; 

	float diffusion_rate <- 0.8;
	int frequence_input <- 3;
	
	int rain_duration <- 60; // Rain stops after N cycles (0 = infinite)
	float infiltration_rate <- 0.001; // Water intake per cycle
	
	list<point> drain_cells <- [];
	list<point> source_cells <- [];
	list<point> points <- flow points_in shape;
	map<point, list<point>> neighbors <- [];
	map<point, float> h <- [];
	float input_water <- 0.02;
	geometry river_g;
	
	list<point> rivers_pt <- [];
	map<point, bool> is_river <- [];
	
	map<point, float> w <- []; // filled elevation
	map<point, list<point>> down <- [];
	map<point, list<float>> down_w <- [];
	float eps <- 0.001;
	
	map<point, bool> is_drain <- [];
	list<point> inf_cells <- [];
	
	float lateral <- 0.02; // sideways conductance on flat cells
	float flat_tol <- 0.01; // elevation tolerance to count as flat
	
	float flood_threshold <- 0.05; // flow height counted as "flooded"
	int flooded_count <- 0;
	
	list<float> level_bins <- [0.0, 0.02, 0.05, 0.1, 0.2, 0.5]; // lower bound of each level
	list<int> level_counts <- [0, 0, 0, 0, 0, 0];
	
	list<float> level_sums <- [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
	int stat_cycles <- 0;	
	list<float> level_avg <- [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

	init {
		flow[] <- 0.0;
		float cell_size <- shape.width / terrain.columns;
		
		// Import street geometry with a 2.0 cell buffer to ensure connected diagonal flow
		create street_import from: file("../includes/gdl_manz_diff.shp");
		if (length(street_import) > 0) {
			river_g <- (union(street_import collect each.shape)) buffer (cell_size * 2.0);
			ask street_import { do die; }
		} else {
			river_g <- shape;
		}
		
		is_river <- points as_map (each::(river_g intersects each));
		rivers_pt <- points where (is_river[each]);
		
		// neighbors from rounded grid indices, street cells only
		float cw <- shape.width / terrain.columns;
		float chh <- shape.height / terrain.rows;
		float x0 <- min(points collect each.x);
		float y0 <- min(points collect each.y);
		int stride <- terrain.columns + 4; // margin so keys never wrap between rows
		map<int, point> by_idx <- [];
		map<point, int> cidx <- [];
		loop p over: rivers_pt {
			int k <- int(round((p.x - x0) / cw)) + int(round((p.y - y0) / chh)) * stride;
			by_idx[k] <- p;
			cidx[p] <- k;
		}
		write "street cells " + length(rivers_pt) + " unique idx " + length(by_idx);
		loop p over: rivers_pt {
			int k <- cidx[p];
			list<point> ln <- [];
			loop d over: [-1, 1, -stride, stride, -stride - 1, -stride + 1, stride - 1, stride + 1] {
				if (by_idx contains_key (k + d)) { add by_idx[k + d] to: ln; }
			}
			neighbors[p] <- ln;
		}
		write "avg neighbors " + mean(neighbors.values collect length(each));
		
		list<float> valid_elevs <- (points collect float(terrain[each])) where (each > 0.0);
		float min_val <- length(valid_elevs) > 0 ? min(valid_elevs) : 0.0;
		
		if (use_dem_elevation) {
			h <- points as_map (each::(float(terrain[each]) <= 0.0 ? min_val : float(terrain[each])));
		} else {
			h <- points as_map (each::0.0);
		}
		
		if (fill) {
			loop pt over: rivers_pt {
				flow[pt] <- 0.2;
			}
		}
		
		if (length(rivers_pt) > 0) {
			float max_h <- max(rivers_pt collect h[each]);
			float min_h <- min(rivers_pt collect h[each]);
			
			// Lowest 6% elevation street cells act as map outlets/sinks
			drain_cells <- rivers_pt where (h[each] <= min_h + ((max_h - min_h) * 0.06));
			
			loop d over: drain_cells { is_drain[d] <- true; }
			do fill_from(drain_cells);
			
			// disconnected street pieces: one outlet at their lowest cell, single pass
			loop s over: rivers_pt {
				if (w contains_key s) { continue; }
				list<point> comp <- [s];
				map<point, bool> seen <- [s::true];
				int i <- 0;
				loop while: i < length(comp) {
					loop n over: neighbors[comp[i]] where (is_river[each] and !(seen contains_key each)) {
						seen[n] <- true;
						add n to: comp;
					}
					i <- i + 1;
				}
				point o <- comp with_min_of h[each];
				add o to: drain_cells;
				is_drain[o] <- true;
				do fill_from([o]);
			}
			
			inf_cells <- rivers_pt where (!(is_drain contains_key each));
			if (rainfall_mode) {
				source_cells <- inf_cells;
			} else {
				source_cells <- rivers_pt where (h[each] >= max_h - ((max_h - min_h) * 0.10));
			}
			
			// static downhill neighbors and slope weights
			loop p over: inf_cells {
				list<point> ln <- neighbors[p] where ((w contains_key each) and w[each] <= w[p] + flat_tol);
				if (!empty(ln)) {
					list<float> s <- ln collect (max(w[p] - w[each], 0.0) / (p distance_to each) + lateral);
					float t <- sum(s);
					down[p] <- ln;
					down_w[p] <- s collect (each / t);
				}
			}
			write "streets " + length(rivers_pt) + " drains " + length(drain_cells) + " down " + length(down);
		}
	}

	// Rain injection onto street network
	reflex adding_input_water when: every(frequence_input#cycle) and (rain_duration <= 0 or cycle <= rain_duration) {
		loop p over: source_cells {
			flow[p] <- flow[p] + input_water;
		}
	}

	// Dynamic Downhill Flow Physics (3 Sub-passes per cycle for fast visible movement)
	reflex flowing {
		loop pass from: 1 to: 3 {
			map<point, float> q <- (down.keys where (flow[each] > 0.0001)) as_map (each::flow[each]);
				loop p over: q.keys {
				float moved <- q[p] * diffusion_rate;
				if (moved > 0.0001) {
					list<point> ln <- down[p];
					list<float> ws <- down_w[p];
					flow[p] <- flow[p] - moved;
					loop i from: 0 to: length(ln) - 1 {
						flow[ln[i]] <- flow[ln[i]] + moved * ws[i];
					}
				}
			}
		}
	}

	// Ground absorption / storm drain intake
	reflex infiltration when: infiltration_rate > 0.0 {
		loop p over: inf_cells {
			if (flow[p] > 0.0) {
				flow[p] <- max(0.0, flow[p] - infiltration_rate);
			}
		}
	}

	// Instant outlet discharge at lowest map boundaries
	reflex draining {
		loop p over: drain_cells {
			flow[p] <- 0.0;
		}
	}
	
	reflex flood_stats {
		flooded_count <- length(rivers_pt where (flow[each] >= flood_threshold));
	}
	
	reflex level_stats {
		list<int> counts <- [];
		loop i from: 0 to: length(level_bins) - 1 {
			float lo <- level_bins[i];
			float hi <- (i = length(level_bins) - 1) ? 99999.0 : level_bins[i + 1];
			add length(rivers_pt where (flow[each] >= lo and flow[each] < hi)) to: counts;
		}
		level_counts <- counts;
		stat_cycles <- stat_cycles + 1;
		loop i from: 0 to: length(level_bins) - 1 {
			level_sums[i] <- level_sums[i] + counts[i];
		}
		loop i from: 0 to: length(level_bins) - 1 {
			level_avg[i] <- level_sums[i] / max(1, stat_cycles);
		}
	}
	
	action fill_from (list<point> seeds) {
		list<point> frontier <- copy(seeds);
		loop s over: seeds { w[s] <- h[s]; }
		loop while: !empty(frontier) {
			list<point> next <- remove_duplicates(frontier accumulate (neighbors[each] where (is_river[each] and !(w contains_key each))));
			loop n over: next {
				w[n] <- max(h[n], min((neighbors[n] where (w contains_key each)) collect w[each]) + eps);
			}
			frontier <- next;
		}
	}

}

species street_import {
}

experiment hydro type: gui {
	parameter "Input water rate" var: input_water <- 0.02 min: 0.001 max: 0.5 step: 0.005;
	parameter "Rainfall Mode (All Streets)" var: rainfall_mode <- true;
	parameter "Rain Duration (Cycles, 0 = continuous)" var: rain_duration <- 60 min: 0 max: 500 step: 10;
	parameter "Infiltration Rate" var: infiltration_rate <- 0.001 min: 0.0 max: 0.05 step: 0.001;
	parameter "Fill streets initially" var: fill <- false;
	parameter "Use DEM Elevation" var: use_dem_elevation <- true;
	parameter "Flood threshold" var: flood_threshold <- 0.05 min: 0.001 max: 0.5 step: 0.005;

	output {
		display d type: 3d {
			mesh terrain scale: 0 triangulation: true color: palette([#burlywood, #saddlebrown, #darkgreen, #green]) refresh: false smooth: true;
			
			graphics "Street Vector Layer" {
				draw river_g color: #yellow wireframe: true;
			}
			
			mesh flow scale: 1.0 triangulation: true color: palette(reverse(brewer_colors("Blues"))) transparency: 0.3 no_data: 0.0;
		}
		display charts refresh: every(1#cycle) {
			chart "Flooded cells" type: series {
				data "flooded" value: flooded_count color: #red;
			}
		}
		display levels refresh: every(1#cycle) {
			chart "Cells by water level" type: histogram {
				data "0.00-0.02" value: level_counts[0] color: #lightblue;
				data "0.02-0.05" value: level_counts[1] color: #deepskyblue;
				data "0.05-0.10" value: level_counts[2] color: #dodgerblue;
				data "0.10-0.20" value: level_counts[3] color: #royalblue;
				data "0.20-0.50" value: level_counts[4] color: #mediumblue;
				data "0.50+"     value: level_counts[5] color: #darkblue;
			}
		}
		display average refresh: every(1#cycle) {
			chart "Average cells by level (all cycles)" type: pie {
				data "0.00-0.02" value: level_avg[0] color: #lightblue;
				data "0.02-0.05" value: level_avg[1] color: #deepskyblue;
				data "0.05-0.10" value: level_avg[2] color: #dodgerblue;
				data "0.10-0.20" value: level_avg[3] color: #royalblue;
				data "0.20-0.50" value: level_avg[4] color: #mediumblue;
				data "0.50+"     value: level_avg[5] color: #darkblue;
			}
		}
				
	}
}