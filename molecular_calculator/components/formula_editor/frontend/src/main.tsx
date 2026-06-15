import React from "react";
import { createRoot } from "react-dom/client";
import FormulaEditor from "./FormulaEditor";
import "./index.css";

createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <FormulaEditor />
  </React.StrictMode>
);
