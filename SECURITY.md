# Security Policy

## Supported Versions

We release patches for security vulnerabilities for the following versions:

| Version | Supported          |
| ------- | ------------------ |
| 0.2.x   | :white_check_mark: |
| < 0.2   | :x:                |

## Reporting a Vulnerability

We take security vulnerabilities seriously. If you discover a security issue, please follow these steps:

### 1. Do Not Open a Public Issue

Please **do not** open a public GitHub issue for security vulnerabilities, as this could put users at risk.

### 2. Report Privately

Send your report privately to the Reactome team:

- **Email**: help@reactome.org
- **Subject**: [SECURITY] Logic Network Generator - Brief description

### 3. Include in Your Report

Please include as much information as possible:

- **Type of vulnerability** (e.g., SQL injection, command injection, XSS)
- **Full paths of affected source files**
- **Location of the affected code** (tag/branch/commit or direct URL)
- **Step-by-step instructions to reproduce** the issue
- **Proof of concept or exploit code** (if possible)
- **Impact of the vulnerability** (what an attacker could do)
- **Suggested fix** (if you have one)

### 4. What to Expect

- **Acknowledgment**: We'll acknowledge receipt of your report within 48 hours
- **Assessment**: We'll assess the vulnerability and determine severity
- **Timeline**: We'll provide an expected timeline for a fix
- **Updates**: We'll keep you informed of progress
- **Credit**: If you wish, we'll credit you in the security advisory

### 5. Disclosure Policy

- We'll work with you to understand and resolve the issue
- We'll aim to patch critical vulnerabilities within 30 days
- We'll coordinate disclosure timing with you
- We'll publicly disclose once a patch is available

## Security Best Practices for Users

### Environment Variables

- Never commit `.env` files or credentials to version control
- Use `.env.example` as a template (never put real credentials here)
- Keep Neo4j connection strings secure

### Neo4j Database

- Use authentication for Neo4j in production
- Don't expose Neo4j ports publicly
- Keep Neo4j version up to date
- Use Docker network isolation when running in containers

### Dependencies

- Regularly update dependencies: `poetry update`
- Check for known vulnerabilities: `poetry show --outdated`
- Review security advisories for dependencies

### Input Validation

- Validate pathway IDs before processing
- Be cautious with pathway lists from untrusted sources
- Sanitize file paths to prevent directory traversal

### Generated Files

- Be careful when sharing generated network files
- They may contain sensitive biological data
- Follow your organization's data handling policies

## Known Security Considerations

### 1. Neo4j Connection

The tool connects to a Neo4j database. Ensure:
- Database connection uses authentication
- Connection string is stored securely (environment variables, not code)
- Database is not publicly accessible

### 2. Command Injection

The tool uses subprocess calls for git operations. We:
- Sanitize all inputs
- Use parameterized commands
- Avoid shell=True where possible

### 3. File System Access

The tool reads from and writes to the file system. Users should:
- Run with minimal necessary permissions
- Restrict output directory permissions
- Validate file paths from external sources

### 4. Dependency Vulnerabilities

We monitor dependencies for known vulnerabilities:
- All dependencies are managed through Poetry
- We use GitHub Dependabot for automated updates
- Security advisories are reviewed promptly

## Vulnerability Disclosure

When a vulnerability is fixed, we will:

1. Release a patch version
2. Publish a GitHub Security Advisory
3. Update CHANGELOG.md with security fix notes
4. Credit the reporter (if they wish)
5. Notify users through release notes

## Security Update Process

1. **Assessment**: Verify and assess the vulnerability
2. **Fix Development**: Develop and test the fix
3. **Testing**: Ensure fix works and doesn't break functionality
4. **Release**: Create a patch release
5. **Notification**: Notify users via GitHub release
6. **Documentation**: Update security documentation

## Contact

For security-related questions or concerns:

- **Email**: help@reactome.org
- **GitHub**: https://github.com/reactome/logic-network-generator/security

## Attribution

This security policy is based on best practices from:
- [GitHub Security Best Practices](https://docs.github.com/en/code-security)
- [OWASP Security Guidelines](https://owasp.org/)
