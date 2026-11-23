# ChromaCode 项目约定

本文件面向所有 AI 编码工具，自包含。`CLAUDE.md` 通过 `@AGENTS.md` 导入本文件，不重复正文。

## 这是什么仓库

面向研究与对外说明的**文档仓库**：`README.md`（英文）、`README.zh.md`（中文）、`assets/`（图与演示）。代码与实现不在本仓库，在 [`ChromaCode-core`](../ChromaCode-core)。

因此本仓库的改动几乎都是文字与图片：没有构建、没有测试，重点是内容准确。

## 工作区布局（tower）

| 位置 | 分支 | 职责 |
| --- | --- | --- |
| `/workspaces/projects/ChromaCode` | `main` | 主仓 |
| `/workspaces/projects/_worktrees/ChromaCode/<任务名>` | `<type>/<任务名>` | 需要多步改写时的工作树 |

改一两句话直接在 `main` 上改并提交即可；只有需要并行推进多个改动时才开 worktree。worktree 不放仓库内部。

## 内容规则

- **中英文同步**：改 `README.md` 必须同步 `README.zh.md`（反之亦然），术语一致。
- **不夸大结论**：区分"仿真结果"与"实测精度"、"研究方向"与"已验证部署"。仓库里已有的 `Research status` 段落体现这一要求，改文字不要削弱它。
- 图片放 `assets/en/`、`assets/zh/` 对应语言目录，命名用小写连字符。
- 不写个人机器路径、账号、私有环境细节。

## 提交

完成并验证一处内容改动后，及时提交并推送。一个提交只做一件事；提交信息写清改了什么、为什么。不擅自提交无关的用户改动；提交或推送受阻时说明原因。不对 `main` 强推。
