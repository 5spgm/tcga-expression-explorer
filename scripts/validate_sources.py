#!/usr/bin/env python3
"""docs/sources.js の構造検証。

sources.js は「文字列と入れ子のオブジェクト/配列だけ」で書かれたデータ
モジュールなので、JS処理系を持ち込まずに検証できる。行うのは3つ:

  1. 文字列リテラルの外側で括弧が釣り合っているか
  2. 閉じ忘れた文字列リテラルがないか
  3. SOURCE_DEFS / DATA_SOURCES / SUBTYPE_SOURCES を実際に組み立てて、
     - 各 source 定義に ja と en の両方があるか
     - DATA_SOURCES の各行が参照する source が定義済みか

使い方:  python3 validate_sources.py docs/sources.js
終了コード 0 で合格、1 で不合格。
"""
import ast
import re
import sys


def strip_strings(src):
    """文字列リテラルを空白に置き換えた版と、リテラルの位置を返す。"""
    out, i, n = [], 0, len(src)
    spans = []
    while i < n:
        ch = src[i]
        if ch in "\"'`":
            q, start, i = ch, i, i + 1
            while i < n and src[i] != q:
                i += 2 if src[i] == "\\" else 1
            if i >= n:
                raise SyntaxError(f"閉じていない文字列リテラル (offset {start})")
            spans.append((start, i))
            out.append(" " * (i - start + 1))
            i += 1
        elif src.startswith("//", i):
            j = src.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i)); i = j
        elif src.startswith("/*", i):
            j = src.find("*/", i + 2)
            if j < 0:
                raise SyntaxError("閉じていないブロックコメント")
            out.append(" " * (j + 2 - i)); i = j + 2
        else:
            out.append(ch); i += 1
    return "".join(out), spans


def check_brackets(bare):
    pairs = {")": "(", "]": "[", "}": "{"}
    stack = []
    for pos, ch in enumerate(bare):
        if ch in "([{":
            stack.append((ch, pos))
        elif ch in pairs:
            if not stack or stack[-1][0] != pairs[ch]:
                raise SyntaxError(f"括弧の不一致: {ch!r} (offset {pos})")
            stack.pop()
    if stack:
        ch, pos = stack[-1]
        raise SyntaxError(f"閉じていない {ch!r} (offset {pos})")


def strip_comments(src):
    """文字列リテラルを壊さずにコメントだけ落とす。

    URL中の '//' を行コメントと誤認しないよう、文字列の内外を追跡する。
    """
    out, i, n = [], 0, len(src)
    while i < n:
        ch = src[i]
        if ch in "\"'`":
            q, start, i = ch, i, i + 1
            while i < n and src[i] != q:
                i += 2 if src[i] == "\\" else 1
            out.append(src[start:i + 1]); i += 1
        elif src.startswith("//", i):
            j = src.find("\n", i)
            i = n if j < 0 else j
        elif src.startswith("/*", i):
            j = src.find("*/", i + 2)
            i = n if j < 0 else j + 2
        else:
            out.append(ch); i += 1
    return "".join(out)


def split_strings(src):
    """(コード片, 文字列リテラル) の並びに分解する。文字列は末尾に来る。"""
    parts, i, n = [], 0, len(src)
    buf = []
    while i < n:
        ch = src[i]
        if ch in "\"'`":
            q, start, i = ch, i, i + 1
            while i < n and src[i] != q:
                i += 2 if src[i] == "\\" else 1
            parts.append(("".join(buf), src[start:i + 1]))
            buf = []
            i += 1
        else:
            buf.append(ch); i += 1
    if buf:
        parts.append(("".join(buf), ""))
    return parts


def js_to_python(block):
    """`key: "a" + "b"` 形式のJSリテラルをPythonの式に直す。

    裸のキーの引用符付けは、文字列リテラルの外側にだけ適用する
    (URLや本文中のコロンを壊さないため)。
    """
    block = strip_comments(block)
    key = re.compile(r'(?<=[{,\[\s])([A-Za-z_$][\w$]*)\s*:')
    out = []
    for code, lit in split_strings(block):
        out.append(key.sub(r'"\1":', code))
        out.append(lit)
    return "".join(out)


_JS_CONST = {"true": True, "false": False, "null": None, "undefined": None}


def safe_eval(node):
    """文字列・数値・真偽値・配列・オブジェクト・文字列連結のみを評価する。"""
    if isinstance(node, ast.Expression):
        return safe_eval(node.body)
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        if node.id in _JS_CONST:
            return _JS_CONST[node.id]
        raise ValueError(f"未対応の識別子: {node.id}")
    if isinstance(node, ast.List):
        return [safe_eval(e) for e in node.elts]
    if isinstance(node, ast.Dict):
        return {safe_eval(k): safe_eval(v) for k, v in zip(node.keys, node.values)}
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left, right = safe_eval(node.left), safe_eval(node.right)
        if isinstance(left, str) or isinstance(right, str):
            return f"{left}{right}"
        return left + right
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        v = safe_eval(node.operand)
        return -v if isinstance(node.op, ast.USub) else v
    raise ValueError(f"未対応のノード: {type(node).__name__}")


def load(src):
    """SOURCE_DEFS / DATA_SOURCES / SUBTYPE_SOURCES を組み立てて返す。"""
    got = {}
    for name in ("SOURCE_DEFS", "DATA_SOURCES", "SUBTYPE_SOURCES"):
        m = re.search(rf"\b(?:const|let|var)\s+{name}\s*=\s*", src)
        if not m:
            continue
        start = m.end()
        bare, _ = strip_strings(src[start:])
        depth, end = 0, None
        for pos, ch in enumerate(bare):
            if ch in "{[":
                depth += 1
            elif ch in "}]":
                depth -= 1
                if depth == 0:
                    end = pos + 1
                    break
        if end is None:
            raise SyntaxError(f"{name} のリテラルが閉じていない")
        py = js_to_python(src[start:start + end])
        got[name] = safe_eval(ast.parse(py, mode="eval"))
    return got


def main(path):
    src = open(path, encoding="utf-8").read()
    problems = []

    bare, spans = strip_strings(src)
    check_brackets(bare)
    print(f"構文OK  文字列リテラル {len(spans)} 個 / {len(src)} 文字")

    ns = load(src)
    sd = ns.get("SOURCE_DEFS")
    if sd is None:
        problems.append("SOURCE_DEFS が見つからない")
    else:
        print(f"SOURCE_DEFS: {len(sd)} 定義")
        for key, v in sd.items():
            for lang in ("ja", "en"):
                if not v.get(lang):
                    problems.append(f"{key}: {lang} が無い")
            if not v.get("name"):
                problems.append(f"{key}: name が無い")

    ds = ns.get("DATA_SOURCES")
    if ds is not None and sd is not None:
        n_rows = 0
        for cancer, rows in ds.items():
            for row in rows:
                n_rows += 1
                s = row.get("source")
                if s and s not in sd:
                    problems.append(f"DATA_SOURCES[{cancer}]: 未定義の source {s!r}")
        print(f"DATA_SOURCES: {len(ds)} がん種 / {n_rows} 行")

    ss = ns.get("SUBTYPE_SOURCES")
    if ss is not None and sd is not None:
        # 配列でもがん種キーの辞書でも受ける
        if isinstance(ss, dict):
            entries = [(k, v) for k, v in ss.items()]
        else:
            entries = [(row.get("cancer", i) if isinstance(row, dict) else i, row)
                       for i, row in enumerate(ss)]
        for label, row in entries:
            rows = row if isinstance(row, list) else [row]
            for r in rows:
                s = r.get("source") if isinstance(r, dict) else None
                if s and s not in sd:
                    problems.append(f"SUBTYPE_SOURCES[{label}]: 未定義の source {s!r}")
        print(f"SUBTYPE_SOURCES: {len(entries)} 件")

    if problems:
        print("\n不合格:")
        for p in problems:
            print("  -", p)
        return 1
    print("\n合格")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("使い方: validate_sources.py <sources.js>")
    sys.exit(main(sys.argv[1]))
