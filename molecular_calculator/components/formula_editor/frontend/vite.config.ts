import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// base: "./" => relative asset paths so Streamlit can serve the built files
// from the component's static directory regardless of mount path.
export default defineConfig({
  base: "./",
  plugins: [react()],
  build: {
    outDir: "build",
    emptyOutDir: true,
  },
});
