% move the robot azimuth if needed
if current_az ~= az
	disp(sprintf('...setting new azimuth %d', az));
	roboaz = 0.001 * 0.5 * interp1(azs, azdata, az);
	set(roboth.AzimuthCtrl, 'Value', roboaz);
	set(roboth.az_value, 'String', sprintf('%.3f', roboaz));
	robot_tdtazmove(servo, roboaz);
	current_az = az;
end

% move the robot elevation if needed
if current_el ~= el
	disp(sprintf('...setting new elevation %d', el));
	roboel = 0.001 * 0.5 * interp1(els, eldata, el);
	set(roboth.ElevationCtrl, 'Value', roboel);
	set(roboth.el_value, 'String', sprintf('%.3f', roboel));
	robot_tdtelmove(servo, roboel);
	current_el = el;
end

% ask user to move robot to current speaker
disp(sprintf('Move speaker to azimuth %d, elevation %d', az, el));
while ~query_user('Continue')	end
