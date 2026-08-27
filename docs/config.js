// ============================================================
// 環境ごとの設定
// ============================================================
// このファイルだけは、配信環境に合わせて各自が書き換える。
// app.js / i18n.js はそのまま差し替えてよい設計にしてあるので、
// 更新のたびに DATA_BASE_URL を書き戻す必要がない。
//
// 末尾にスラッシュを付けないこと
// (app.js 側が `${DATA_BASE_URL}/manifest.json` の形で組み立てる)。
//
// 配信元には以下が必要:
//   - HTTPS であること(GitHub PagesがHTTPSのため、混在コンテンツになる)
//   - Access-Control-Allow-Origin でサイトのオリジンを許可していること
//
// 例:
//   Cloudflare R2   "https://pub-xxxxxxxx.r2.dev"
//   自前のnginx     "https://example.jp/data"

const DATA_BASE_URL = "https://pub-19864563e9014e228cefc601d77adfbc.r2.dev";
