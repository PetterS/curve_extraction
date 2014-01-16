% Make sure all settings are given on correct format
function settings = parse_settings(settings)

	int_entries = {'num_threads', 'maxiter'};
	for i = 1:numel(int_entries)
		if (isfield(settings,int_entries{i}))
				val = int32(getfield(settings,  int_entries{i}));
				settings = setfield(settings, int_entries{i},val);
		end
	end

	bool_entries = {'verbose','use_a_star','store_visit_time','store_parents' };
	for i = 1:numel(bool_entries)
		if (isfield(settings,bool_entries{i}))
				val = logical(getfield(settings,  bool_entries{i}));
			settings = setfield(settings, bool_entries{i},val);
		end
	end

end