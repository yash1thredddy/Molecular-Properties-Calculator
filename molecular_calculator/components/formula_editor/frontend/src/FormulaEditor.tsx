import React, { ReactNode } from "react";
import {
  Streamlit,
  StreamlitComponentBase,
  withStreamlitConnection,
} from "streamlit-component-lib";
import Editor, { Monaco, loader } from "@monaco-editor/react";

// Pin Monaco core to the version in package.json so the runtime (CDN) version
// can't drift from the build. (Finding #15.) For fully offline/air-gapped use,
// replace this URL with a locally-bundled vs path. (Finding #11.)
loader.config({
  paths: { vs: "https://cdn.jsdelivr.net/npm/monaco-editor@0.52.2/min/vs" },
});

/**
 * Monaco-based formula editor for Streamlit.
 * Args from Python: value (string), columns (string[]), functions (string[]),
 *   height (number).
 * Returns: the current formula string (via Streamlit.setComponentValue).
 *
 * v1 (this spike): local column + function autocomplete via Monaco's
 * registerCompletionItemProvider. v2 will add inline AI completions.
 */
class FormulaEditor extends StreamlitComponentBase {
  private completionDisposable: { dispose: () => void } | null = null;

  public render = (): ReactNode => {
    const value = (this.props.args["value"] as string) ?? "";
    const height = (this.props.args["height"] as number) ?? 160;
    const dark = this.props.theme?.base === "dark";

    return (
      <div className="fe-wrap">
        <Editor
          height={height}
          defaultLanguage="plaintext"
          defaultValue={value}
          theme={dark ? "vs-dark" : "light"}
          options={{
            minimap: { enabled: false },
            lineNumbers: "off",
            fontSize: 14,
            scrollBeyondLastLine: false,
            wordWrap: "on",
            quickSuggestions: true,
            suggestOnTriggerCharacters: true,
            renderLineHighlight: "none",
            folding: false,
            overviewRulerLanes: 0,
            // Don't auto-insert a closing "]" when the user types "[" — the
            // column completion already inserts the full "[Column]", and
            // auto-closing would produce "[[Column]]".
            autoClosingBrackets: "never",
            autoSurround: "never",
          }}
          onMount={this.handleMount}
          onChange={this.handleChange}
        />
      </div>
    );
  };

  public componentDidUpdate = (): void => {
    // Re-register completions if columns/functions changed across reruns.
    if (this.monaco) this.registerCompletions(this.monaco);
    Streamlit.setFrameHeight();
  };

  private monaco: Monaco | null = null;

  private handleMount = (_editor: unknown, monaco: Monaco): void => {
    this.monaco = monaco;
    this.registerCompletions(monaco);
    Streamlit.setFrameHeight();
  };

  private registerCompletions = (monaco: Monaco): void => {
    if (this.completionDisposable) this.completionDisposable.dispose();
    const columns = (this.props.args["columns"] as string[]) ?? [];
    const functions = (this.props.args["functions"] as string[]) ?? [];

    this.completionDisposable = monaco.languages.registerCompletionItemProvider(
      "plaintext",
      {
        triggerCharacters: ["[", "(", " "],
        provideCompletionItems: (model: any, position: any) => {
          const word = model.getWordUntilPosition(position);
          const lineContent = model.getLineContent(position.lineNumber);
          let startColumn = word.startColumn;
          let endColumn = word.endColumn;
          // If the user already typed "[" just before the word, consume it so the
          // inserted "[Column]" doesn't double into "[[Column]]". (Columns are
          // 1-indexed; lineContent is 0-indexed, so column c is lineContent[c-1].)
          if (startColumn > 1 && lineContent[startColumn - 2] === "[") {
            startColumn -= 1;
          }
          // Consume a trailing "]" sitting at the cursor for the same reason.
          if (lineContent[endColumn - 1] === "]") {
            endColumn += 1;
          }
          const range = {
            startLineNumber: position.lineNumber,
            endLineNumber: position.lineNumber,
            startColumn,
            endColumn,
          };
          const cols = columns.map((c) => ({
            label: `[${c}]`,
            kind: monaco.languages.CompletionItemKind.Field,
            insertText: `[${c}]`,
            detail: "column",
            range,
          }));
          const fns = functions.map((f) => ({
            label: `${f}()`,
            kind: monaco.languages.CompletionItemKind.Function,
            insertText: `${f}($0)`,
            insertTextRules:
              monaco.languages.CompletionItemInsertTextRule.InsertAsSnippet,
            detail: "function",
            range,
          }));
          return { suggestions: [...cols, ...fns] };
        },
      }
    );
  };

  private handleChange = (value: string | undefined): void => {
    Streamlit.setComponentValue(value ?? "");
  };
}

export default withStreamlitConnection(FormulaEditor);
