rsqrd = mkSyntheticImage([65 65],'rsqrd');
displayImage(rsqrd);

gaussian = mkSyntheticImage([65 65],'gaussian',[-4 4],[-4 4],2);
displayImage(gaussian);

disc=mkDisc([64 64], 16);
displayImage(disc);

fractal=mkFract([64 64], 2);
displayImage(fractal);

gauss=mkGaussian([64 64]);
displayImage(gauss);

impulse=mkImpulse([64 64]);
displayImage(impulse);

rad=mkR([64 64]);
displayImage(rad);

ramp=mkRamp([64 64],pi/4);
displayImage(ramp);

sine=mkSine([64 64],16,0);
displayImage(sine);

square=mkSquare([64 64],8,0);
displayImage(square);

theta=mkAngle([64 64]);
displayImage(theta);

zone=mkZonePlate([64 64]);
displayImage(zone);

