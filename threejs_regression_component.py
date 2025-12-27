"""
Three.js 3D OLS Regression Viewer Component

This module provides an interactive 3D visualization for OLS regression using Three.js.
It offers superior interactivity compared to Plotly's 3D charts with features like:
- Smooth rotation, zoom, and pan
- Orthographic/perspective camera toggle
- Show/hide residual vectors
- Show/hide predicted markers on plane
- Color points by residual magnitude
- Interactive axis arrows
- Click on points to trigger callbacks

Developed by: Yashwanth Reddy for ITR-UIC
Part of: Chemo-Informatics Toolkit
"""

import json
import numpy as np
from typing import Optional, List, Dict, Any


def get_threejs_regression_component(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    z_predicted: np.ndarray,
    residuals: np.ndarray,
    b0: float,
    b1: float,
    b2: float,
    r_squared: float,
    x_label: str = "X",
    y_label: str = "Y",
    z_label: str = "Z",
    chart_id: str = "threejs_regression",
    smiles_data: Optional[List[str]] = None,
    name_data: Optional[List[str]] = None,
    height: int = 700
) -> str:
    """
    Generate Three.js 3D OLS regression visualization component.

    This component provides an interactive 3D visualization with:
    - Data points colored by residual magnitude
    - Fitted plane surface
    - Optional residual vectors
    - Optional predicted markers
    - Axis arrows with labels
    - Camera controls (orbit, zoom, pan)
    - Perspective/orthographic toggle

    Args:
        x: X values (independent variable 1)
        y: Y values (independent variable 2)
        z: Z values (dependent variable - actual)
        z_predicted: Predicted Z values on the fitted plane
        residuals: Residual values (z - z_predicted)
        b0: Intercept coefficient
        b1: X coefficient
        b2: Y coefficient
        r_squared: R-squared value
        x_label: Label for X axis
        y_label: Label for Y axis
        z_label: Label for Z axis
        chart_id: Unique identifier for this chart
        smiles_data: Optional list of SMILES strings for each point
        name_data: Optional list of names/IDs for each point
        height: Height of the component in pixels

    Returns:
        HTML string containing the complete Three.js component
    """
    # Convert numpy arrays to lists for JSON serialization
    x_list = x.tolist() if isinstance(x, np.ndarray) else list(x)
    y_list = y.tolist() if isinstance(y, np.ndarray) else list(y)
    z_list = z.tolist() if isinstance(z, np.ndarray) else list(z)
    z_pred_list = z_predicted.tolist() if isinstance(z_predicted, np.ndarray) else list(z_predicted)
    residuals_list = residuals.tolist() if isinstance(residuals, np.ndarray) else list(residuals)

    # Prepare SMILES and name data
    smiles_json = json.dumps(smiles_data if smiles_data else [])
    names_json = json.dumps(name_data if name_data else [])

    # Calculate data ranges for scaling
    x_min, x_max = float(np.min(x)), float(np.max(x))
    y_min, y_max = float(np.min(y)), float(np.max(y))
    z_min, z_max = float(np.min(z)), float(np.max(z))

    # Add padding
    x_pad = 0.1 * (x_max - x_min) if x_max > x_min else 1.0
    y_pad = 0.1 * (y_max - y_min) if y_max > y_min else 1.0
    z_pad = 0.1 * (z_max - z_min) if z_max > z_min else 1.0

    # Equation string
    equation = f"{z_label} = {b0:.3f} + {b1:.3f}*{x_label} + {b2:.3f}*{y_label}"

    # This component runs entirely within the Streamlit iframe
    html_component = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            overflow: hidden;
        }}
        #container {{
            width: 100%;
            height: {height}px;
            position: relative;
            background: linear-gradient(135deg, #0f0f1a 0%, #1a1a2e 50%, #16213e 100%);
            border-radius: 8px;
            overflow: hidden;
        }}
        #canvas3d {{
            width: 100%;
            height: 100%;
            display: block;
        }}
        .controls-panel {{
            position: absolute;
            top: 10px;
            left: 10px;
            background: rgba(255,255,255,0.95);
            padding: 12px 16px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            font-size: 13px;
            z-index: 100;
            max-width: 220px;
        }}
        .controls-panel h4 {{
            margin: 0 0 10px 0;
            padding-bottom: 8px;
            border-bottom: 1px solid #eee;
            color: #333;
            font-size: 14px;
        }}
        .control-item {{
            display: flex;
            align-items: center;
            margin-bottom: 8px;
            cursor: pointer;
        }}
        .control-item input {{
            margin-right: 8px;
        }}
        .control-item label {{
            cursor: pointer;
            user-select: none;
        }}
        .reset-btn {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 12px;
            width: 100%;
            margin-top: 10px;
        }}
        .reset-btn:hover {{
            opacity: 0.9;
        }}
        .stats-panel {{
            position: absolute;
            top: 10px;
            right: 10px;
            background: rgba(255,255,255,0.95);
            padding: 12px 16px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            font-family: 'Monaco', 'Menlo', monospace;
            font-size: 12px;
            z-index: 100;
            max-width: 320px;
        }}
        .stats-panel h4 {{
            margin: 0 0 8px 0;
            color: #333;
            font-weight: 600;
        }}
        .stats-panel .r2 {{
            color: #667eea;
            font-weight: 600;
        }}
        .stats-panel .equation {{
            color: #666;
            font-size: 11px;
            word-break: break-all;
            margin-top: 4px;
        }}
        .stats-panel .count {{
            color: #999;
            font-size: 10px;
            margin-top: 8px;
        }}
        .tooltip {{
            position: absolute;
            display: none;
            background: rgba(0,0,0,0.9);
            color: white;
            padding: 10px 14px;
            border-radius: 6px;
            font-size: 12px;
            pointer-events: none;
            z-index: 1000;
            max-width: 280px;
            border: 1px solid rgba(255,255,255,0.2);
        }}
        .instructions {{
            position: absolute;
            bottom: 10px;
            left: 50%;
            transform: translateX(-50%);
            background: rgba(0,0,0,0.6);
            color: white;
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 11px;
            z-index: 100;
        }}
        .loading {{
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            color: white;
            font-size: 16px;
            z-index: 50;
        }}
        .axis-legend {{
            position: absolute;
            bottom: 50px;
            right: 10px;
            background: rgba(0,0,0,0.7);
            padding: 8px 12px;
            border-radius: 6px;
            font-size: 11px;
            color: white;
            z-index: 100;
        }}
        .axis-legend div {{
            margin: 3px 0;
            display: flex;
            align-items: center;
            gap: 6px;
        }}
        .axis-legend .color-box {{
            width: 12px;
            height: 12px;
            border-radius: 2px;
        }}
    </style>
</head>
<body>
    <div id="container">
        <canvas id="canvas3d"></canvas>

        <div class="loading" id="loading">Loading Three.js...</div>

        <div class="controls-panel" id="controls" style="display: none;">
            <h4>3D Visualization Controls</h4>

            <div class="control-item">
                <input type="checkbox" id="orthoToggle">
                <label for="orthoToggle">Orthographic projection</label>
            </div>

            <div class="control-item">
                <input type="checkbox" id="axisToggle" checked>
                <label for="axisToggle">Show axis arrows</label>
            </div>

            <div class="control-item">
                <input type="checkbox" id="gridToggle" checked>
                <label for="gridToggle">Show grid lines</label>
            </div>

            <div class="control-item">
                <input type="checkbox" id="residualToggle">
                <label for="residualToggle">Show residual vectors</label>
            </div>

            <div class="control-item">
                <input type="checkbox" id="predictedToggle">
                <label for="predictedToggle">Show predicted markers</label>
            </div>

            <div class="control-item">
                <input type="checkbox" id="autoRotateToggle">
                <label for="autoRotateToggle">Auto-rotate</label>
            </div>

            <button class="reset-btn" id="resetBtn">Reset Camera</button>
        </div>

        <div class="stats-panel" id="stats" style="display: none;">
            <h4>Model Statistics</h4>
            <div>R² = <span class="r2">{r_squared:.4f}</span></div>
            <div class="equation">{equation}</div>
            <div class="count">n = {len(x_list)} points</div>
        </div>

        <div class="tooltip" id="tooltip"></div>

        <div class="axis-legend" id="axisLegend" style="display: none;">
            <div><span class="color-box" style="background:#ff4444;"></span> X: {x_label}</div>
            <div><span class="color-box" style="background:#4488ff;"></span> Z: {z_label}</div>
            <div><span class="color-box" style="background:#44ff44;"></span> Y: {y_label}</div>
        </div>

        <div class="instructions">
            Drag to rotate | Scroll to zoom | Right-click to pan | Hover for details
        </div>
    </div>

    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>

    <script>
    (function() {{
        'use strict';

        // Data
        const xData = {json.dumps(x_list)};
        const yData = {json.dumps(y_list)};
        const zData = {json.dumps(z_list)};
        const zPredicted = {json.dumps(z_pred_list)};
        const residualsData = {json.dumps(residuals_list)};
        const smilesData = {smiles_json};
        const namesData = {names_json};

        // Original data ranges (without padding for tick labels)
        const xMinOrig = {x_min}, xMaxOrig = {x_max};
        const yMinOrig = {y_min}, yMaxOrig = {y_max};
        const zMinOrig = {z_min}, zMaxOrig = {z_max};

        // Data ranges with padding
        const xMin = {x_min - x_pad}, xMax = {x_max + x_pad};
        const yMin = {y_min - y_pad}, yMax = {y_max + y_pad};
        const zMin = {z_min - z_pad}, zMax = {z_max + z_pad};

        // Plane coefficients
        const b0 = {b0}, b1 = {b1}, b2 = {b2};

        // Labels
        const xLabel = "{x_label}";
        const yLabel = "{y_label}";
        const zLabel = "{z_label}";

        // Wait for Three.js to load
        function waitForThree() {{
            if (typeof THREE !== 'undefined' && typeof THREE.OrbitControls !== 'undefined') {{
                document.getElementById('loading').style.display = 'none';
                document.getElementById('controls').style.display = 'block';
                document.getElementById('stats').style.display = 'block';
                document.getElementById('axisLegend').style.display = 'block';
                initScene();
            }} else {{
                setTimeout(waitForThree, 100);
            }}
        }}

        waitForThree();

        function initScene() {{
            const container = document.getElementById('container');
            const canvas = document.getElementById('canvas3d');
            const tooltip = document.getElementById('tooltip');

            const width = container.clientWidth;
            const height = container.clientHeight;

            // Normalize function to map data to -5 to 5 range
            function normalize(val, min, max) {{
                const range = max - min || 1;
                return (val - min) / range * 10 - 5;
            }}

            // Scene - slightly lighter background for better label visibility
            const scene = new THREE.Scene();
            scene.background = new THREE.Color(0x1a1a2e);

            // Camera - positioned for better view of axes
            let camera = new THREE.PerspectiveCamera(55, width / height, 0.1, 1000);
            camera.position.set(15, 12, 15);
            camera.lookAt(0, 0, 0);

            // Renderer
            const renderer = new THREE.WebGLRenderer({{
                canvas: canvas,
                antialias: true
            }});
            renderer.setSize(width, height);
            renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

            // Controls
            const controls = new THREE.OrbitControls(camera, renderer.domElement);
            controls.enableDamping = true;
            controls.dampingFactor = 0.05;
            controls.autoRotate = false;
            controls.autoRotateSpeed = 1.0;
            controls.target.set(0, 0, 0);

            // Lighting
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.7);
            scene.add(ambientLight);

            const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
            directionalLight.position.set(10, 20, 10);
            scene.add(directionalLight);

            const directionalLight2 = new THREE.DirectionalLight(0xffffff, 0.4);
            directionalLight2.position.set(-10, -10, -10);
            scene.add(directionalLight2);

            // Groups for toggling visibility
            const axisGroup = new THREE.Group();
            const gridGroup = new THREE.Group();
            const residualGroup = new THREE.Group();
            const predictedGroup = new THREE.Group();

            scene.add(axisGroup);
            scene.add(gridGroup);
            scene.add(residualGroup);
            scene.add(predictedGroup);

            // Store point meshes for raycasting
            const pointMeshes = [];

            // Color scale for residuals
            const maxResidual = Math.max(...residualsData.map(Math.abs));

            function getResidualColor(residual) {{
                const normalized = Math.abs(residual) / (maxResidual || 1);
                const r = Math.min(1, 0.3 + normalized * 0.7);
                const g = Math.max(0, 0.85 - normalized * 0.75);
                const b = Math.max(0, 0.4 - normalized * 0.3);
                return new THREE.Color(r, g, b);
            }}

            // Create data points
            const pointGeometry = new THREE.SphereGeometry(0.2, 20, 20);

            for (let i = 0; i < xData.length; i++) {{
                const material = new THREE.MeshPhongMaterial({{
                    color: getResidualColor(residualsData[i]),
                    shininess: 80
                }});

                const sphere = new THREE.Mesh(pointGeometry, material);
                sphere.position.set(
                    normalize(xData[i], xMin, xMax),
                    normalize(zData[i], zMin, zMax),  // Y is up in Three.js
                    normalize(yData[i], yMin, yMax)
                );

                sphere.userData = {{
                    index: i,
                    x: xData[i],
                    y: yData[i],
                    z: zData[i],
                    zPred: zPredicted[i],
                    residual: residualsData[i],
                    smiles: smilesData[i] || null,
                    name: namesData[i] || null
                }};

                scene.add(sphere);
                pointMeshes.push(sphere);
            }}

            // Create fitted plane using PlaneGeometry
            const planeSize = 10;
            const planeSegments = 20;
            const planeGeometry = new THREE.PlaneGeometry(planeSize, planeSize, planeSegments, planeSegments);

            // Modify vertices to match the fitted plane equation
            const positions = planeGeometry.attributes.position;
            for (let i = 0; i < positions.count; i++) {{
                const px = positions.getX(i);  // -5 to 5
                const pz = positions.getY(i);  // -5 to 5 (becomes Z in world)

                // Convert normalized coords back to data space
                const dataX = (px + 5) / 10 * (xMax - xMin) + xMin;
                const dataY = (pz + 5) / 10 * (yMax - yMin) + yMin;

                // Calculate Z using plane equation: z = b0 + b1*x + b2*y
                const dataZ = b0 + b1 * dataX + b2 * dataY;
                const normalizedZ = normalize(dataZ, zMin, zMax);

                // Set the Y position (up axis) to the calculated Z
                positions.setZ(i, normalizedZ);
            }}

            planeGeometry.computeVertexNormals();

            const planeMaterial = new THREE.MeshPhongMaterial({{
                color: 0x4a90d9,
                transparent: true,
                opacity: 0.55,
                side: THREE.DoubleSide,
                shininess: 30
            }});

            const plane = new THREE.Mesh(planeGeometry, planeMaterial);
            plane.rotation.x = -Math.PI / 2;
            scene.add(plane);

            // ============= SIMPLIFIED 3D AXES WITH BRIGHT LABELS =============

            // Helper function to create text sprite with background for visibility
            function createTextSprite(text, fontSize, color, withBackground) {{
                const canvas2d = document.createElement('canvas');
                const ctx = canvas2d.getContext('2d');
                canvas2d.width = 512;
                canvas2d.height = 256;

                // Add semi-transparent dark background for better readability
                if (withBackground) {{
                    ctx.fillStyle = 'rgba(0, 0, 0, 0.7)';
                    const padding = 20;
                    ctx.beginPath();
                    ctx.roundRect(padding, padding, canvas2d.width - 2*padding, canvas2d.height - 2*padding, 15);
                    ctx.fill();
                }}

                // Draw text with bright color and outline for visibility
                ctx.font = `bold ${{fontSize || 64}}px Arial, sans-serif`;
                ctx.textAlign = 'center';
                ctx.textBaseline = 'middle';

                // Add text outline/stroke for better visibility
                ctx.strokeStyle = '#000000';
                ctx.lineWidth = 6;
                ctx.strokeText(text, 256, 128);

                // Fill with bright color
                ctx.fillStyle = color || '#ffffff';
                ctx.fillText(text, 256, 128);

                const texture = new THREE.CanvasTexture(canvas2d);
                texture.minFilter = THREE.LinearFilter;
                const spriteMaterial = new THREE.SpriteMaterial({{
                    map: texture,
                    transparent: true,
                    depthTest: false  // Always render on top
                }});
                const sprite = new THREE.Sprite(spriteMaterial);
                sprite.scale.set(3, 1.5, 1);
                return sprite;
            }}

            // Helper function to format tick values
            function formatTickValue(val) {{
                if (Math.abs(val) >= 1000) return val.toExponential(1);
                if (Math.abs(val) >= 100) return val.toFixed(0);
                if (Math.abs(val) >= 10) return val.toFixed(1);
                if (Math.abs(val) >= 1) return val.toFixed(2);
                return val.toFixed(3);
            }}

            // Create axis line with tick marks and BRIGHT labels
            function createAxis(startPos, endPos, color, label, dataMin, dataMax, numTicks) {{
                numTicks = numTicks || 5;

                // Convert hex color to CSS string for sprites
                const colorHex = '#' + color.toString(16).padStart(6, '0');

                // Main axis line (thicker)
                const axisMaterial = new THREE.LineBasicMaterial({{ color: color, linewidth: 3 }});
                const axisGeometry = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(...startPos),
                    new THREE.Vector3(...endPos)
                ]);
                const axisLine = new THREE.Line(axisGeometry, axisMaterial);
                axisGroup.add(axisLine);

                // Calculate axis direction
                const dir = new THREE.Vector3(
                    endPos[0] - startPos[0],
                    endPos[1] - startPos[1],
                    endPos[2] - startPos[2]
                );
                const length = dir.length();
                dir.normalize();

                // Create tick marks and labels
                const tickMaterial = new THREE.LineBasicMaterial({{ color: 0x888888 }});

                for (let i = 0; i <= numTicks; i++) {{
                    const t = i / numTicks;
                    const pos = new THREE.Vector3(
                        startPos[0] + dir.x * length * t,
                        startPos[1] + dir.y * length * t,
                        startPos[2] + dir.z * length * t
                    );

                    // Tick mark (perpendicular line)
                    const tickLength = 0.3;
                    let perpDir;
                    if (Math.abs(dir.y) > 0.9) {{
                        // Vertical axis - tick goes outward in X-Z plane
                        perpDir = new THREE.Vector3(-1, 0, -1).normalize();
                    }} else if (Math.abs(dir.x) > 0.5) {{
                        // X axis - tick goes down and back
                        perpDir = new THREE.Vector3(0, -1, -1).normalize();
                    }} else {{
                        // Z axis (our Y data) - tick goes down and left
                        perpDir = new THREE.Vector3(-1, -1, 0).normalize();
                    }}

                    const tickStart = pos.clone();
                    const tickEnd = pos.clone().add(perpDir.clone().multiplyScalar(tickLength));

                    const tickGeometry = new THREE.BufferGeometry().setFromPoints([tickStart, tickEnd]);
                    const tick = new THREE.Line(tickGeometry, tickMaterial);
                    axisGroup.add(tick);

                    // Tick label with actual data value - BRIGHT WHITE with background
                    const dataValue = dataMin + (dataMax - dataMin) * t;
                    const tickLabel = createTextSprite(formatTickValue(dataValue), 48, '#FFFFFF', false);
                    tickLabel.scale.set(1.2, 0.6, 1);
                    tickLabel.position.copy(tickEnd).add(perpDir.clone().multiplyScalar(0.6));
                    axisGroup.add(tickLabel);
                }}

                // Axis label at the end - BRIGHT COLOR with background
                const labelSprite = createTextSprite(label, 72, colorHex, true);
                labelSprite.scale.set(3.5, 1.75, 1);
                const labelPos = new THREE.Vector3(...endPos).add(dir.clone().multiplyScalar(2.0));
                labelSprite.position.copy(labelPos);
                axisGroup.add(labelSprite);
            }}

            // Create the three main axes with tick marks - BRIGHT COLORS
            // X axis (bright red) - from left to right along X
            createAxis([-5, -5, -5], [5, -5, -5], 0xff4444, xLabel, xMinOrig, xMaxOrig, 5);

            // Z axis in Three.js (bright green, represents Y data) - from front to back
            createAxis([-5, -5, -5], [-5, -5, 5], 0x44ff44, yLabel, yMinOrig, yMaxOrig, 5);

            // Y axis in Three.js (bright blue, represents Z data) - vertical
            createAxis([-5, -5, -5], [-5, 5, -5], 0x4488ff, zLabel, zMinOrig, zMaxOrig, 5);

            // ============= 3 SIDE GRIDS ONLY (like scipy plot) =============

            const gridMaterial = new THREE.LineBasicMaterial({{ color: 0x555577, transparent: true, opacity: 0.4 }});
            const gridDivisions = 8;

            // Floor grid (XY plane at z=-5, which is Y in Three.js terms)
            for (let i = 0; i <= gridDivisions; i++) {{
                const t = (i / gridDivisions) * 10 - 5;

                // Lines parallel to X
                const geom1 = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(-5, -5, t),
                    new THREE.Vector3(5, -5, t)
                ]);
                gridGroup.add(new THREE.Line(geom1, gridMaterial));

                // Lines parallel to Z (Y in data)
                const geom2 = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(t, -5, -5),
                    new THREE.Vector3(t, -5, 5)
                ]);
                gridGroup.add(new THREE.Line(geom2, gridMaterial));
            }}

            // Back wall grid (at z=-5 in Three.js, XZ plane in data coords)
            for (let i = 0; i <= gridDivisions; i++) {{
                const t = (i / gridDivisions) * 10 - 5;

                // Horizontal lines
                const geom1 = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(-5, t, -5),
                    new THREE.Vector3(5, t, -5)
                ]);
                gridGroup.add(new THREE.Line(geom1, gridMaterial));

                // Vertical lines
                const geom2 = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(t, -5, -5),
                    new THREE.Vector3(t, 5, -5)
                ]);
                gridGroup.add(new THREE.Line(geom2, gridMaterial));
            }}

            // Left wall grid (at x=-5, YZ plane in data coords)
            for (let i = 0; i <= gridDivisions; i++) {{
                const t = (i / gridDivisions) * 10 - 5;

                // Horizontal lines
                const geom1 = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(-5, t, -5),
                    new THREE.Vector3(-5, t, 5)
                ]);
                gridGroup.add(new THREE.Line(geom1, gridMaterial));

                // Depth lines
                const geom2 = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(-5, -5, t),
                    new THREE.Vector3(-5, 5, t)
                ]);
                gridGroup.add(new THREE.Line(geom2, gridMaterial));
            }}

            // Create residual vectors (initially hidden)
            residualGroup.visible = false;
            const lineMaterial = new THREE.LineBasicMaterial({{ color: 0xf1c40f, linewidth: 2 }});

            for (let i = 0; i < xData.length; i++) {{
                const points = [
                    new THREE.Vector3(
                        normalize(xData[i], xMin, xMax),
                        normalize(zData[i], zMin, zMax),
                        normalize(yData[i], yMin, yMax)
                    ),
                    new THREE.Vector3(
                        normalize(xData[i], xMin, xMax),
                        normalize(zPredicted[i], zMin, zMax),
                        normalize(yData[i], yMin, yMax)
                    )
                ];
                const geometry = new THREE.BufferGeometry().setFromPoints(points);
                const line = new THREE.Line(geometry, lineMaterial);
                residualGroup.add(line);
            }}

            // Create predicted markers (initially hidden)
            predictedGroup.visible = false;
            const predGeometry = new THREE.SphereGeometry(0.1, 8, 8);
            const predMaterial = new THREE.MeshBasicMaterial({{ color: 0x19d3f3 }});

            for (let i = 0; i < xData.length; i++) {{
                const marker = new THREE.Mesh(predGeometry, predMaterial);
                marker.position.set(
                    normalize(xData[i], xMin, xMax),
                    normalize(zPredicted[i], zMin, zMax),
                    normalize(yData[i], yMin, yMax)
                );
                predictedGroup.add(marker);
            }}

            // Raycaster for hover detection
            const raycaster = new THREE.Raycaster();
            const mouse = new THREE.Vector2();
            let hoveredPoint = null;

            function onMouseMove(event) {{
                const rect = canvas.getBoundingClientRect();
                mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
                mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

                raycaster.setFromCamera(mouse, camera);
                const intersects = raycaster.intersectObjects(pointMeshes);

                if (intersects.length > 0) {{
                    const point = intersects[0].object.userData;
                    canvas.style.cursor = 'pointer';

                    // Highlight point
                    if (hoveredPoint !== intersects[0].object) {{
                        if (hoveredPoint) {{
                            hoveredPoint.material.emissive.setHex(0x000000);
                            hoveredPoint.scale.set(1, 1, 1);
                        }}
                        hoveredPoint = intersects[0].object;
                        hoveredPoint.material.emissive.setHex(0x333333);
                        hoveredPoint.scale.set(1.3, 1.3, 1.3);
                    }}

                    let tooltipHTML = `<div style="margin-bottom: 6px; font-weight: 600; color: #8be9fd;">Point ${{point.index + 1}}</div>`;

                    if (point.name) {{
                        tooltipHTML += `<div><strong>Name:</strong> ${{point.name}}</div>`;
                    }}

                    tooltipHTML += `
                        <div><strong>${{xLabel}}:</strong> ${{point.x.toFixed(3)}}</div>
                        <div><strong>${{yLabel}}:</strong> ${{point.y.toFixed(3)}}</div>
                        <div><strong>${{zLabel}}:</strong> ${{point.z.toFixed(3)}}</div>
                        <div><strong>Predicted:</strong> ${{point.zPred.toFixed(3)}}</div>
                        <div style="color: #ff6b6b;"><strong>Residual:</strong> ${{point.residual.toFixed(3)}}</div>
                    `;

                    tooltip.innerHTML = tooltipHTML;
                    tooltip.style.display = 'block';
                    tooltip.style.left = (event.clientX - rect.left + 15) + 'px';
                    tooltip.style.top = (event.clientY - rect.top + 15) + 'px';
                }} else {{
                    canvas.style.cursor = 'grab';
                    tooltip.style.display = 'none';

                    if (hoveredPoint) {{
                        hoveredPoint.material.emissive.setHex(0x000000);
                        hoveredPoint.scale.set(1, 1, 1);
                        hoveredPoint = null;
                    }}
                }}
            }}

            canvas.addEventListener('mousemove', onMouseMove);
            canvas.addEventListener('mouseleave', () => {{
                tooltip.style.display = 'none';
                if (hoveredPoint) {{
                    hoveredPoint.material.emissive.setHex(0x000000);
                    hoveredPoint.scale.set(1, 1, 1);
                    hoveredPoint = null;
                }}
            }});

            // Control handlers
            document.getElementById('orthoToggle').addEventListener('change', function() {{
                const aspect = width / height;
                const frustumSize = 15;

                if (this.checked) {{
                    const newCam = new THREE.OrthographicCamera(
                        frustumSize * aspect / -2,
                        frustumSize * aspect / 2,
                        frustumSize / 2,
                        frustumSize / -2,
                        0.1, 1000
                    );
                    newCam.position.copy(camera.position);
                    newCam.rotation.copy(camera.rotation);
                    camera = newCam;
                }} else {{
                    const newCam = new THREE.PerspectiveCamera(60, aspect, 0.1, 1000);
                    newCam.position.copy(camera.position);
                    newCam.rotation.copy(camera.rotation);
                    camera = newCam;
                }}
                controls.object = camera;
                controls.update();
            }});

            document.getElementById('axisToggle').addEventListener('change', function() {{
                axisGroup.visible = this.checked;
            }});

            document.getElementById('gridToggle').addEventListener('change', function() {{
                gridGroup.visible = this.checked;
            }});

            document.getElementById('residualToggle').addEventListener('change', function() {{
                residualGroup.visible = this.checked;
            }});

            document.getElementById('predictedToggle').addEventListener('change', function() {{
                predictedGroup.visible = this.checked;
            }});

            document.getElementById('autoRotateToggle').addEventListener('change', function() {{
                controls.autoRotate = this.checked;
            }});

            document.getElementById('resetBtn').addEventListener('click', function() {{
                camera.position.set(15, 12, 15);
                camera.lookAt(0, 0, 0);
                controls.target.set(0, 0, 0);
                controls.update();
            }});

            // Animation loop
            function animate() {{
                requestAnimationFrame(animate);
                controls.update();
                renderer.render(scene, camera);
            }}

            animate();

            // Handle resize
            window.addEventListener('resize', function() {{
                const newWidth = container.clientWidth;
                const newHeight = container.clientHeight;
                camera.aspect = newWidth / newHeight;
                camera.updateProjectionMatrix();
                renderer.setSize(newWidth, newHeight);
            }});

            console.log('[Three.js Regression] Scene initialized with', xData.length, 'points');
        }}
    }})();
    </script>
</body>
</html>
    """

    return html_component


def get_threejs_colorbar() -> str:
    """
    Generate HTML for a color bar legend showing residual magnitude scale.

    Returns:
        HTML string for the color bar
    """
    return """
    <div style="
        display: flex;
        align-items: center;
        gap: 10px;
        margin-top: 10px;
        padding: 10px;
        background: #f8f9fa;
        border-radius: 6px;
    ">
        <span style="font-size: 12px; color: #666;">|Residual|:</span>
        <div style="
            width: 150px;
            height: 15px;
            background: linear-gradient(to right, #4dcc4d, #ffcc00, #ff4444);
            border-radius: 3px;
        "></div>
        <span style="font-size: 11px; color: #999;">Low -> High</span>
    </div>
    """
