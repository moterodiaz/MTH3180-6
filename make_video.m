function make_video()
cases = string_modes();
writer = VideoWriter('vibrating_string.mp4','MPEG-4');
writer.FrameRate = 30;
writer.Quality = 100;
open(writer);
cleanup = onCleanup(@() close(writer));
for k = 1:numel(cases)
    s = cases(k);
    fprintf('Rendering n=%d, mode %d, omega=%.6f rad/s\n',s.p.n,s.mode,s.omega);
    animate_string(s.t,s.V,s.p,'',writer,s.mode,s.omega);
end
end
